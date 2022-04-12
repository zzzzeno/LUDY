#!/usr/bin/env python3

import sys
import os
import numpy as np
pjoin = os.path.join

import multiprocessing
import observable as obs
import vegas
import pickle
from NLO_integrands import NLO_integrands
import time

from pprint import pprint, pformat

import vegas
import argparse

import logging
logger = logging.getLogger("LU_DY")
logger.setLevel(logging.INFO)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)

class MockUpRes(object):
    def __init__(self, res, jac=1., j1=None, j2=None, pg=None, spin1=None, spin2=None):
        self.res = res
        self.jac = jac
        self.j1 = j1
        self.j2 = j2
        self.pg = pg
        self.spin1 = spin1
        self.spin2 = spin2

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

class LU_DY(object):

    def __init__(self, *args, 
                    mZ=91.188, eCM=13000., resolution_coll=0., resolution_soft=0., h_sigma=2, 
                    initial_state_a=0., min_pt=10., max_pt=2000., n_bins=1, basis=None, tag="tqq", mode='min_pt', phase='real',
                    verbosity = 0, n_cores=None
                ):

        if basis is None:
            # Use a default basis
            basis = np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))

        if n_cores is None:
            self.n_cores = multiprocessing.cpu_count()
        else:
            self.n_cores = n_cores

        self.verbosity = verbosity

        self.integrand = NLO_integrands(
            mZ, eCM, resolution_coll, resolution_soft, h_sigma, 
            initial_state_a, min_pt, max_pt, n_bins, basis, 
            tag=tag, debug=(self.verbosity>2), mode=mode, phase=phase
        )

        self.tag = tag
        self.observables = None
        self.worker_observables = None
        self.pool = None
        self.timings = None

    # Helper function for paralellisation
    @staticmethod
    def evaluate_sample(args):

        i_job, integrand, observables, timings, apply_observables, call_args = args
        x, integrator_wgt, h_cube = call_args
        
        if i_job is not None:
            return i_job, (observables if apply_observables else None), [ LU_DY.evaluate_sample((None, integrand, observables, timings, apply_observables, (a_x,a_wgt,a_h_cube) )) for a_x, a_wgt, a_h_cube in zip(x, integrator_wgt, h_cube) ]

        start = time.time()
        
        all_res = integrand.safe_eval(x)

        # The mockup function below is useful for validation
        #all_res = [ MockUpRes( sum(xi**2 for xi in x) ), ]

        if timings is not None:
            timings.increment('c++',time.time()-start) 

        final_weight = 0.
        for r in all_res:
            final_weight += r.res*r.jac*(1. if integrator_wgt is None else integrator_wgt)

        evt_group = None
        if observables is not None:
            start = time.time()
            evt_group = obs.EventGroup([
                obs.Event(
                        initial_state_jets = obs.JetList( [ obs.Jet(r.j1, spin=r.spin1), obs.Jet(r.j2, spin=r.spin2) ] ), 
                        final_state_jets = obs.JetList( [ obs.Jet(r.pg, spin=3 ), ] ),
                        weight = r.res*r.jac*(1. if integrator_wgt is None else integrator_wgt),
                        E_com = integrand.eCM
                ) for r in all_res
            ], h_cube=h_cube)
            for e in evt_group:
                e.compute_derived_quantities()
            if apply_observables:
                observables.accumulate_event_group(evt_group)
                evt_group = None
            if timings is not None:
                timings.increment('observables', time.time()-start)

        return (final_weight, evt_group)

    @vegas.batchintegrand
    def evaluate_samples(self, xs, wgts=None, h_cubes=None):    
        
        if self.n_cores == 1:
            if wgts is None:
                return [ LU_DY.evaluate_sample((None, self.integrand, self.worker_observables[0], self.timings, True, (x,None,None) ))[0] for x in xs ]
            else:
                return [ LU_DY.evaluate_sample((None, self.integrand, self.worker_observables[0], self.timings, True, (x,wgt,h_cube) ))[0] for x, wgt, h_cube in zip(xs, wgts, h_cubes) ]
        else:
             
            x_chunks = list(chunks(xs,1+len(xs)//self.n_cores))
            wgts_chunks = list(chunks(wgts,1+len(xs)//self.n_cores)) if wgts is not None else [[None,]*len(x_chunk) for x_chunk in x_chunks]
            h_cubes = list(chunks(h_cubes,1+len(xs)//self.n_cores)) if h_cubes is not None else [[None,]*len(x_chunk) for x_chunk in x_chunks]

            if self.verbosity>0:
                logger.info("Dispatching %d points over %d cores."%(len(xs), self.n_cores))

            # It's more efficient to apply the observables in parallel since they are now no longer acting on shared memory so there is no locking anymore.
            apply_observable_in_parallel = True
            jobs_it = self.pool.imap(LU_DY.evaluate_sample, 
                    list( (i_job, self.integrand, self.worker_observables[i_job], self.timings, apply_observable_in_parallel, job_args) for i_job, job_args in enumerate(zip(x_chunks, wgts_chunks, h_cubes)) ))

            all_res = {}
            for i_job, observables, res in jobs_it:
                all_res[i_job] = [r[0] for r in res]
                if observables is not None:
                    self.worker_observables[i_job] = observables
                if self.verbosity>0:
                    if not apply_observable_in_parallel:
                        print("N jobs done: %d/%d, now processing events..."%(len(all_res),len(x_chunks)), end='\r')
                    else:
                        print("N jobs done: %d/%d"%(len(all_res),len(x_chunks)), end='\r')
                if not apply_observable_in_parallel:
                    start = time.time()
                    for r in res:
                        self.worker_observables[i_job].accumulate_event_group(r[1])
                    timings.increment('observables', time.time()-start)
                    print("N jobs done: %d/%d, events processed."%(len(all_res),len(x_chunks)), end='\r')
            
            return sum([all_res[k] for k in sorted(list(all_res.keys()))],[])

    def integrate(self, n_iterations=10, n_evals_training=10000, n_evals_production=10000, vegas_grid='./vegas_grid.pkl', 
            alpha=0.5, beta=0.75, nhcube_batch=1000, seed=None, target_relative_error=0., target_absolute_error=0., 
            observables=None, hwu_path=None, **opts
        ):
        
        if seed:
            np.random.seed(seed)
        
        if hwu_path is None:
            hwu_path = 'histograms_%s.HwU'%self.tag
        else:
            if not hwu_path.endswith('HwU'):
                hwu_path = '%s.HwU'%hwu_path

        starting_map = None
        if vegas_grid is not None and vegas_grid.upper()!='NONE' and os.path.isfile(vegas_grid):
            logger.info("Recycling trained VEGAS grid from '%s'."%vegas_grid)
            starting_map = pickle.load(open(vegas_grid,'rb'))

        integrator = vegas.Integrator(
            map=([(0.,1.),]*9+[(0.,1.),] if starting_map is None else starting_map),
            nitn=(n_iterations if starting_map is None else 1),
            neval=(n_evals_training if starting_map is None else n_evals_production),
            alpha=alpha,
            beta=beta,
            adapt=(starting_map is None),
            analyzer=(vegas.reporter() if (self.verbosity>0) else None),
            nhcube_batch=nhcube_batch,
            rtol=target_relative_error,
            atol=target_absolute_error,
            mpi=False,
            **opts
        )

        if starting_map is None:
            # Disable observables during training
            self.observables = None
            self.worker_observables = [None for _ in range(self.n_cores)]
            integration_start = time.time()
            logger.info("Training VEGAS grid with %d iterations of %d points..."%(n_iterations, n_evals_training))
            with obs.LUDYManager() as manager:
                with multiprocessing.Pool(processes=self.n_cores) as self.pool:
                    this_manager = (manager if self.n_cores>1 else None)
                    if this_manager is None:
                        self.timings = obs.Timings()
                    else:
                        self.timings = manager.Timings()
                    training_res = integrator(self.evaluate_samples, nitn=n_iterations, neval=n_evals_training, adapt=True)
                    dict_timings = self.timings.to_dict()
                    tot_CPU = dict_timings['c++']+dict_timings['observables']
                    wall_time = time.time()-integration_start
                    logger.info("Result after training:\n%s"%str(training_res.summary()))
                    logger.info("Integration wall time %.0f s. tot CPU time %.0f (x%.1f) (C++ : %.3g%%, Python observables: %.3g%% )"%(
                        wall_time,
                        tot_CPU,
                        tot_CPU/wall_time,
                        (dict_timings['c++']/tot_CPU)*100.,
                        (dict_timings['observables']/tot_CPU)*100.,
                    ))
            if vegas_grid is not None and vegas_grid.upper()!='NONE':
                logger.info("Saving trained VEGAS grid to '%s'."%vegas_grid)
                pickle.dump(integrator.map,open(vegas_grid,'wb'))

        if n_evals_production>0:

            # Now build observables
            with obs.LUDYManager() as manager:
                with multiprocessing.Pool(processes=self.n_cores) as self.pool:

                    integration_start = time.time()

                    this_manager = (manager if self.n_cores>1 else None)

                    if this_manager is None:
                        self.timings = obs.Timings()
                    else:
                        self.timings = manager.Timings()

                    # Enable observables
                    self.observables = obs.ObservableList([ obs.CrossSection(), ])
                    self.worker_observables = [obs.ObservableList([ obs.CrossSection(), ]) for _ in range(self.n_cores)]
                    if observables is not None:
                        for obs_string in observables:
                            if obs_string == 'x1':
                                self.observables.append( obs.x1() )
                                for w_obs in self.worker_observables:
                                    w_obs.append( obs.x1() )
                            else:
                                raise NotImplementedError("Observable %s not implemented yet."%obs_string)

                    logger.info("Starting production run with %d sample points."%n_evals_production)

                    integrator.set(nitn=1, neval=n_evals_production, adapt=False)
                    for xs, wgts, hcubes in integrator.random_batch(yield_hcube=True):
                        if self.verbosity>0:
                            logger.info("Submitting new batch of %d points..."%len(xs))
                        self.evaluate_samples(xs, wgts=wgts, h_cubes=hcubes)

                    for w_obs in self.worker_observables:
                        self.observables.finalize_iteration(worker_observable=w_obs)
                    logger.info("Production results: %.5g +/- %.3g (n_events = %d)"%(
                        self.observables[0].histogram.bins[1].integral,
                        self.observables[0].histogram.bins[1].variance**0.5,
                        self.observables[0].histogram.bins[1].n_entries
                    ))
                    dict_timings = self.timings.to_dict()
                    tot_CPU = dict_timings['c++']+dict_timings['observables']
                    wall_time = time.time()-integration_start
                    logger.info("Integration wall time %.0f s. tot CPU time %.0f (x%.1f) (C++ : %.3g%%, Python observables: %.3g%% )"%(
                        wall_time,
                        tot_CPU,
                        tot_CPU/wall_time,
                        (dict_timings['c++']/tot_CPU)*100.,
                        (dict_timings['observables']/tot_CPU)*100.,
                    ))
                    with open(hwu_path,'w') as f:
                        f.write(self.observables.format_to_HwU())
                    logger.info("Histograms output to file '%s'. They can be rendered using the madgraph/various/histograms.py script."%hwu_path)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Compute NLO Drell-Yan cross-section with Local Unitarity.""")
    requiredNamed = parser.add_argument_group('required named arguments')
    # Required options
    #requiredNamed.add_argument(...)
    # Optional options

    parser.add_argument('--verbosity', '-v', dest='verbosity', type=int, default=1,
                        help='Specify LU DY verbosity (default: %(default)s).')
    parser.add_argument('--seed', '-seed', dest='seed', type=int, default=None,
                        help='Specify random seed (default: %(default)s).')
    parser.add_argument('--batch_size', '-bs', dest='batch_size', type=int, default=1000,
                        help='Batch sizse (default: %(default)s).')
    parser.add_argument('--hwu_path', '-hwu', dest='hwu_path', type=str, default=None,
                        help='Path to output hwu file (default: as per tag).')
    parser.add_argument('--n_evals_training', '-net', dest='n_evals_training', type=int, default=10000,
                        help='Number of sample points for training (default: %(default)s).')
    parser.add_argument('--n_evals_production', '-nep', dest='n_evals_production', type=int, default=10000,
                        help='Number of sample points for production (default: %(default)s).')
    parser.add_argument('--n_iterations_training', '-nit', dest='n_iterations_training', type=int, default=10,
                        help='Number of iterations for training VEGAS grids (default: %(default)s).')
    parser.add_argument('--n_cores', '-c', dest='n_cores', type=int, default=multiprocessing.cpu_count(),
                        help='Specify number of cores for the parallelisation (default: %(default)s).')
    parser.add_argument('--observables', '-obs', dest='observables', type=str, nargs='?', default=None,
                        help='List observables to consider. Use keyword "ALL" to select them all (default: %(default)s).')
    parser.add_argument('--tag', '-tag', dest='tag', type=str, default='tqq', choices=('st','tqg','s','tqq','u'),
                        help='Specify which integrand to run (default: %(default)s).')
    parser.add_argument('--vegas_grid', '-vg', dest='vegas_grid', type=str, default='./vegas_grid.pkl',
                        help='Specify path to save/load VEGAS grid. Specify "none" to disable that and use keyword "clean" to forcefully remove default one at "vegas_grid.pkl". (default: %(default)s).')
    args = parser.parse_args()

    if args.vegas_grid == 'clean':
        if os.path.isfile('./vegas_grid.pkl'):
            logger.info("Removing existing grid 'vegas_grid.pkl' to start a fresh new integration.")
            os.remove('./vegas_grid.pkl')
        args.vegas_grid = './vegas_grid.pkl'

    lu_dy = LU_DY(
        verbosity=args.verbosity,
        n_cores=args.n_cores,
        tag = args.tag
    )

    lu_dy.integrate(
        observables=args.observables,
        seed = args.seed,
        vegas_grid = args.vegas_grid,
        nhcube_batch = args.batch_size,
        hwu_path = args.hwu_path,
        n_evals_training = args.n_evals_training,
        n_iterations = args.n_iterations_training,
        n_evals_production = args.n_evals_production
    )
