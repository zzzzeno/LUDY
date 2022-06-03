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
from integrand_wrapper import set_r, set_kin, set_sig, set_defo, set_MUV
import time
import copy
import math

from pprint import pprint, pformat

import vegas
import argparse

import logging
logger = logging.getLogger("LU_DY")
logger.setLevel(logging.INFO)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)

_CONVERT_TO_FLAV = { 0: 21, -1: -1, -2: -2, 1: 1, 2: 2 }

class MockUpRes(object):
    def __init__(self, res, jac=1., j1=None, j2=None, pg=None, spin1=0, spin2=0):
        self.res = res
        self.jac = jac
        if j1 is None:
            self.j1 = [100.,0.,100.,0.]
        else:
            self.j1 = j1
        if j2 is None:
            self.j2 = [100.,0.,100.,0.]
        else:
            self.j2 = j2
        if pg is None:
            self.pg = [100.,0.,100.,0.]
        else:
            self.pg = pg
        self.spin1 = spin1
        self.spin2 = spin2

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Helper function for parallelising the final combination of observables
def finalize_obs_helper(args):
    job_id, job_args = args
    master_obs, worker_obs = job_args
    res = []
    for m_o, w_o in zip(master_obs, worker_obs):
        m_o.finalize_iteration(worker_observable=obs.ObservableList.merge(w_o))
        res.append(m_o)

    return (job_id, res)

class LU_DY(object):

    def __init__(self, *args, 
                    mZ=91.188, eCM=300., resolution_coll=0., resolution_soft=0., h_sigma=2, 
                    initial_state_a=0., min_pt=10., max_pt=2000., max_x1=1., max_x2=1., n_bins=1, basis=None, tag="tqq", mode='min_pt', phase='real',
                    verbosity = 0, n_cores=None, separate_h_cube_errors=True, digits_f64=8, digits_f128=4, **opts 
                ):

        if basis is None:
            # Use a default basis
            basis = np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))

        if n_cores is None:
            self.n_cores = multiprocessing.cpu_count()
        else:
            self.n_cores = n_cores

        self.verbosity = verbosity
        self.separate_h_cube_errors = separate_h_cube_errors
        self.selector_variables = {'min_pt': min_pt, 'max_pt': max_pt, 'max_x1': max_x1, 'max_x2': max_x2}
        
        self.integrand = [ NLO_integrands(
            mZ, eCM, resolution_coll, resolution_soft, h_sigma, 
            initial_state_a, min_pt, max_pt, n_bins, basis, 
            tag=tag, debug=(self.verbosity>2), mode=mode, phase=phase, obs_tag='JADE'
        ) for _ in range(self.n_cores) ]
        for itg in self.integrand:
            itg.set_digits( digits_f64, digits_f128 )
        set_kin(mZ,eCM)
        set_r(resolution_coll, resolution_soft)
        set_sig(h_sigma)
        # TODO Zeno add all other remaining necessary hyperparameters I may not be aware of, like set_MUV, defo, etc...

        self.tag = tag
        self.observables = None
        self.worker_observables = None
        self.pool = None
        self.timings = None

    @staticmethod
    def build_observables(observables):
        """ Helper function to centralise where observables are defined."""

        obs_list = obs.ObservableList([ obs.CrossSection(), ])
        if observables is not None:
            if len(observables) == 1 and observables[0] == 'paper':
                
                # Start with no observables since the inclusive cross-section will already be added later
                obs_list = obs.ObservableList([])
                
                # Observables summed semi-inclusively over flavours
                obs_list.append( obs.x1FixedFlav(flavour=0, take_abs=True) )
                obs_list.append( obs.x1FixedFlav(flavour=1, take_abs=True) )
                obs_list.append( obs.x1FixedFlav(flavour=2, take_abs=True) )
                obs_list.append( obs.x2FixedFlav(flavour=0, take_abs=True) )
                obs_list.append( obs.x2FixedFlav(flavour=1, take_abs=True) )
                obs_list.append( obs.x2FixedFlav(flavour=2, take_abs=True) )

                # Observables also declined per specific flavour
                for flavour_configs in [
                        None,
                        (-2,-2),(-2,-1),(-2,0),(-2,1),(-2,2),
                        (-1,-2),(-1,-1),(-1,0),(-1,1),(-1,2),
                        (0,-2),(0,-1),(0,0),(0,1),(0,2),
                        (1,-2),(1,-1),(1,0),(1,1),(1,2),
                        (2,-2),(2,-1),(2,0),(2,1),(2,2)                    
                    ]:
                    obs_list.append( obs.CrossSection(flavours=flavour_configs) )
                    obs_list.append( obs.x1(flavours=flavour_configs) ) 
                    obs_list.append( obs.x2(flavours=flavour_configs) )
                    
                    obs_list.append( obs.ptj(min_value=0., max_value=300., n_bins=100, flavours=flavour_configs) )
                    obs_list.append( obs.ptj(title='ptjLogY', min_value=0., max_value=300., n_bins=100, y_axis='log', flavours=flavour_configs) )
                    obs_list.append( obs.ptj(title='ptjZoom', min_value=0., max_value=20., n_bins=100, flavours=flavour_configs) )
                    obs_list.append( obs.log10ptj(min_value=-1., max_value=4., n_bins=100, flavours=flavour_configs) )
                    obs_list.append( obs.log10ptj(title='log10ptjLogY', min_value=-1., max_value=4., n_bins=100, y_axis='log', flavours=flavour_configs) )
                    
                    obs_list.append( obs.z(title='z', min_value=0., max_value=1, n_bins=100, y_axis='lin', flavours=flavour_configs) )
                    obs_list.append( obs.z(title='zLogY', min_value=0., max_value=1, n_bins=100, y_axis='log', flavours=flavour_configs) )
                    obs_list.append( obs.z(title='zZoomX1e2', min_value=0., max_value=1.0e-2, n_bins=100, y_axis='lin', flavours=flavour_configs) )
                    obs_list.append( obs.z(title='zZoomX1e2LogY', min_value=0., max_value=1.0e-2, n_bins=100, y_axis='log', flavours=flavour_configs) )
                    obs_list.append( obs.log10z(title='log10z', min_value=-4., max_value=0., n_bins=100, y_axis='lin', flavours=flavour_configs) )
                    obs_list.append( obs.log10z(title='log10zLogY', min_value=-4., max_value=0., n_bins=100, y_axis='log', flavours=flavour_configs) )
                    
                    obs_list.append( obs.rap(title='rap', min_value=0., max_value=100., n_bins=100, y_axis='lin', flavours=flavour_configs) )
                    obs_list.append( obs.rap(title='rapLogY', min_value=0., max_value=100., n_bins=100, y_axis='log', flavours=flavour_configs) )
                    obs_list.append( obs.rap(title='log10rap', min_value=-3., max_value=3., n_bins=100, y_axis='lin', flavours=flavour_configs) )
                    obs_list.append( obs.rap(title='log10rapLogY', min_value=-3., max_value=3., n_bins=100, y_axis='log', flavours=flavour_configs) )

            else:
                if 'ptj' in observables or 'ALL' in observables:
                    obs_list.append( obs.ptj() )
                    obs_list.append( obs.ptj(title='ptjLogY', min_value=0., max_value=2000., n_bins=100, y_axis='log') )
                    obs_list.append( obs.ptj(title='ptjZoom', min_value=0., max_value=20., n_bins=100) )
                if 'log10ptj' in observables or 'ALL' in observables:
                    obs_list.append( obs.log10ptj() )
                    obs_list.append( obs.log10ptj(title='log10ptjLogY', min_value=-1., max_value=4., n_bins=100, y_axis='log') )
                if 'x1' in observables or 'ALL' in observables:
                    obs_list.append( obs.x1() )
                if 'x2' in observables or 'ALL' in observables:
                    obs_list.append( obs.x2() )
                if 'z' in observables or 'ALL' in observables:
                    obs_list.append( obs.z(title='z', min_value=0., max_value=1, n_bins=100, y_axis='lin') )
                    obs_list.append( obs.z(title='zLogY', min_value=0., max_value=1, n_bins=100, y_axis='log') )
                    obs_list.append( obs.z(title='zZoomX1e2', min_value=0., max_value=1.0e-2, n_bins=100, y_axis='lin') )
                    obs_list.append( obs.z(title='zZoomX1e2LogY', min_value=0., max_value=1.0e-2, n_bins=100, y_axis='log') )
                    obs_list.append( obs.z(title='zZoomX1e4', min_value=0., max_value=1.0e-4, n_bins=100, y_axis='lin') )
                    obs_list.append( obs.z(title='zZoomX1e4LogY', min_value=0., max_value=1.0e-4, n_bins=100, y_axis='log') )
                if 'log10z' in observables or 'ALL' in observables:
                    obs_list.append( obs.log10z(title='log10z', min_value=-6., max_value=0., n_bins=100, y_axis='lin') )
                    obs_list.append( obs.log10z(title='log10zLogY', min_value=-6., max_value=0., n_bins=100, y_axis='log') )
                if ('x1_fixed_flavours' in observables or 'ALL' in observables) and not 'no_x1_fixed_flavours' in observables:
                    obs_list.append( obs.x1FixedFlav(flavour=0, take_abs=True) )
                    obs_list.append( obs.x1FixedFlav(flavour=1, take_abs=True) )
                    obs_list.append( obs.x1FixedFlav(flavour=2, take_abs=True) )
                if ('x2_fixed_flavours' in observables or 'ALL' in observables) and not 'no_x2_fixed_flavours' in observables:
                    obs_list.append( obs.x2FixedFlav(flavour=0, take_abs=True) )
                    obs_list.append( obs.x2FixedFlav(flavour=1, take_abs=True) )
                    obs_list.append( obs.x2FixedFlav(flavour=2, take_abs=True) )
        return obs_list

    @staticmethod
    def cuts(evt, selector_variables):
        """ Implements cuts."""

        if evt.x1 > selector_variables['max_x1']:
            return False
        
        if evt.x2 > selector_variables['max_x2']:
            return False

        if len(evt.final_state_jets)==0:
            if selector_variables['min_pt']>0.:
                return False
        else:
            total_final_state_jet_pt = sum(j.p.pt() for j in evt.final_state_jets)
            if not (selector_variables['min_pt'] < total_final_state_jet_pt < selector_variables['max_pt']):
                return False
        
        return True

    # Helper function for paralellisation
    @staticmethod
    def evaluate_sample(args):
        i_job, integrand, observables, timings, separate_h_cube_errors, apply_observables, selector_variables, call_args = args
        x, integrator_wgt, h_cube = call_args
        try:            
            if i_job is not None:
                return i_job, timings, (observables if apply_observables else None), [ LU_DY.evaluate_sample((None, integrand, observables, timings, separate_h_cube_errors, apply_observables, selector_variables, (a_x,a_wgt,a_h_cube) )) for a_x, a_wgt, a_h_cube in zip(x, integrator_wgt, h_cube) ]

            start = time.time()
            
            #TODO streamline and make hyperparameters
            ap=10.*math.tan(math.pi*(x[9]-1/2.))
            a_jac = 10.*math.pi/pow(math.sin(math.pi*x[9]),2)
            integrand.set_a(ap)
            all_res = integrand.safe_eval(x)
            for r in all_res:
               r.jac *= a_jac
            # The mockup function below is useful for validation
            # all_res = [ MockUpRes( sum(xi**2 for xi in x) ), ]

            if timings is not None:
                timings.increment('c++',time.time()-start) 

            evt_group = None

            start = time.time()
            evt_group = obs.EventGroup([
                obs.Event(
                        initial_state_jets = obs.JetList( [ obs.Jet(r.j1, flavour=_CONVERT_TO_FLAV[int(r.spin1)]), obs.Jet(r.j2, flavour=_CONVERT_TO_FLAV[int(r.spin2)]) ] ), 
                        final_state_jets = obs.JetList( [ obs.Jet(r.pg, flavour=None ), ] ),
                        weight = r.res*r.jac*(1. if integrator_wgt is None else integrator_wgt),
                        E_com = integrand.eCM
                ) for r in all_res if abs(r.res*r.jac*(1. if integrator_wgt is None else integrator_wgt))!=0.0
            ], h_cube=(h_cube if separate_h_cube_errors else 1) )
            
            evt_group.compute_derived_quantities()

            # Filter events to apply selector
            evt_group = obs.EventGroup([ e for e in evt_group if LU_DY.cuts(e, selector_variables)])

            final_weight = 0.
            for evt in evt_group:
                final_weight += evt.weight

            # Useful printout of all events passing by
            #print(evt_group)
            if apply_observables and observables is not None:
                observables.accumulate_event_group(evt_group)
                evt_group = None
            if timings is not None:
                timings.increment('observables', time.time()-start)

            return (final_weight, evt_group)
        except KeyboardInterrupt as e:
            logger.warning("Worker #%d interrupted by user."%i_job)
            pass

    @vegas.batchintegrand
    def evaluate_samples(self, xs, wgts=None, h_cubes=None):    
        
        if self.n_cores == 1:
            if wgts is None:
                return [ LU_DY.evaluate_sample((None, self.integrand[0], self.worker_observables[0], self.timings, self.separate_h_cube_errors, True, self.selector_variables, (x,None,None) ))[0] for x in xs ]
            else:
                return [ LU_DY.evaluate_sample((None, self.integrand[0], self.worker_observables[0], self.timings, self.separate_h_cube_errors, True, self.selector_variables, (x,wgt,h_cube) ))[0] for x, wgt, h_cube in zip(xs, wgts, h_cubes) ]
        else:
             
            x_chunks = list(chunks(xs,1+len(xs)//self.n_cores))
            wgts_chunks = list(chunks(wgts,1+len(xs)//self.n_cores)) if wgts is not None else [[None,]*len(x_chunk) for x_chunk in x_chunks]
            h_cubes = list(chunks(h_cubes,1+len(xs)//self.n_cores)) if h_cubes is not None else [[None,]*len(x_chunk) for x_chunk in x_chunks]

            if self.verbosity>2 or (self.verbosity>0 and wgts is None):
                logger.info("Dispatching %d points over %d cores."%(len(xs), self.n_cores))

            # It's more efficient to apply the observables in parallel since they are now no longer acting on shared memory so there is no locking anymore.
            apply_observable_in_parallel = True
            jobs_it = self.pool.imap(LU_DY.evaluate_sample, 
                    list( (i_job, self.integrand[i_job], self.worker_observables[i_job], obs.Timings(), self.separate_h_cube_errors, apply_observable_in_parallel, self.selector_variables, job_args) for i_job, job_args in enumerate(zip(x_chunks, wgts_chunks, h_cubes)) ))

            all_res = {}
            for i_job, timings, observables, res in jobs_it:
                self.timings.aggregate_timings(timings)
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
                    if self.timings is not None: 
                        self.timings.increment('observables', time.time()-start)
                    print("N jobs done: %d/%d, events processed."%(len(all_res),len(x_chunks)), end='\r')
            
            return sum([all_res[k] for k in sorted(list(all_res.keys()))],[])

    def integrate(self, n_iterations_training=10, n_iterations_production=3, n_evals_training=10000, n_evals_production=10000, vegas_grid='./vegas_grid.pkl', 
            alpha=0.5, beta=0.75, nhcube_batch=1000, seed=None, target_relative_error=0., target_absolute_error=0., 
            observables=None, hwu_path=None, use_checkpoints=True, parallelise_final_combination = True, **opts
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
            nitn=(n_iterations_training if starting_map is None else 1),
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
            logger.info("Training VEGAS grid with %d iterations of %d points..."%(n_iterations_training, n_evals_training))
            with multiprocessing.Pool(processes=self.n_cores) as self.pool:
                
                self.timings = obs.Timings()
                training_res = integrator(self.evaluate_samples, nitn=n_iterations_training, neval=n_evals_training, adapt=True)
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

        if n_evals_production>0 and n_iterations_production>0:

            # Now build observables
            with multiprocessing.Pool(processes=self.n_cores) as self.pool:

                integration_start = time.time()

                self.timings = obs.Timings()

                # Enable observables
                self.observables = LU_DY.build_observables(observables)

                self.worker_observables = [copy.deepcopy(self.observables) for _ in range(self.n_cores)]

                logger.info("Starting production run with %d sample points."%n_evals_production)
                integrator.set(nitn=1, neval=n_evals_production, adapt=False)
                n_production_iteration_performed = 0

                interrupted = False
                n_batch_per_iteration = None
                for i_iteration in range(n_iterations_production):
                    for w_o in self.worker_observables:
                        w_o.reset()
                    if i_iteration>0 and use_checkpoints:
                        if self.verbosity>1: logger.info("Saving current results after iteration #%d to disk..."%i_iteration)
                        start=time.time()
                        #pickle.dump( (self.observables, self.worker_observables), open('./current_results.pkl','wb') )
                        pickle.dump( self.observables, open('./current_results.pkl','wb') )
                        if self.verbosity>1: logger.info("Done. [%.0fs]"%(time.time()-start))

                    for i_batch, (xs, wgts, hcubes) in enumerate(integrator.random_batch(yield_hcube=True)):
                        if i_batch == 0:
                            n_batch_per_iteration = (n_evals_production//len(xs)+1)
                        if self.verbosity>0: logger.info("Iteration #%d/%d : Submitting new batch #%d/%d (%.1f%%) of %d points over %d cores..."%(
                            n_production_iteration_performed+1,n_iterations_production,i_batch+1,n_batch_per_iteration, 
                            ((i_batch+1)/float(n_batch_per_iteration))*100., len(xs), self.n_cores))
                        try:
                            self.evaluate_samples(xs, wgts=wgts, h_cubes=hcubes)
                        except KeyboardInterrupt as e:
                            if not use_checkpoints:
                                logger.warning("Integration aborted by user and checkpoints are disabled.")
                                sys.exit(0)
                            if n_production_iteration_performed==0.:
                                logger.warning("Integration aborted by user and no iteration has completed yet.")
                                sys.exit(0)

                            interrupted = True
                            logger.warning("Integration aborted by user. Will proceed with existing data up to iteration #%d."%i_iteration)
                            if self.verbosity>1: logger.info("Loading current results from disk...")
                            start=time.time()
                            #self.observables, self.worker_observables = pickle.load(open('./current_results.pkl','rb'))
                            self.observables = pickle.load(open('./current_results.pkl','rb'))
                            if self.verbosity>1: logger.info("Done. [%.0fs]"%(time.time()-start))
                            break
                        if self.verbosity>1: logger.info("Iteration #%d : Batch evaluation completed."%(n_production_iteration_performed+1))
                    if interrupted:
                        break
                    
                    n_production_iteration_performed += 1
                    if not parallelise_final_combination or self.n_cores==1:
                        logger.info("Finalising iteration #%d/%d and combining %d histograms..."%(n_production_iteration_performed,n_iterations_production,len(self.observables)))
                        self.observables.finalize_iteration(worker_observable=obs.ObservableList.merge(self.worker_observables))
                    else:
                        logger.info("Finalising iteration #%d/%d and combining %d histograms over %d cores..."%(n_production_iteration_performed,n_iterations_production,len(self.observables), self.n_cores))

                        jobs_input_master_obs = [ obs.ObservableList([o,]) for o in self.observables ]
                        jobs_input_worker_obs = [ [ obs.ObservableList([wobs[i_obs],]) for wobs in self.worker_observables ] for i_obs in range(len(self.observables)) ]
                        jobs_input_master_obs_chunks =  list(chunks(jobs_input_master_obs,1+len(jobs_input_master_obs)//self.n_cores))
                        jobs_input_worker_obs_chunks =  list(chunks(jobs_input_worker_obs,1+len(jobs_input_worker_obs)//self.n_cores))

                        jobs_it = self.pool.imap(finalize_obs_helper, list( enumerate(zip(jobs_input_master_obs_chunks, jobs_input_worker_obs_chunks)) ) )
                        all_res = {}
                        for i_job, processes_observable in jobs_it:
                            all_res[i_job] = processes_observable
                            print("N jobs done: %d/%d"%(len(all_res),len(jobs_input_master_obs_chunks)), end='\r')
                        
                        self.observables = obs.ObservableList([ ol[0] for ol in sum([all_res[i] for i in sorted(list(all_res.keys()))],[]) ])

                    tmp_obs = copy.deepcopy(self.observables)
                    tmp_obs.clean_up(normalisation_factor=(1./n_production_iteration_performed))
                    hwu_path_for_this_iteration = '%s_iteration_%d.%s'%('.'.join(hwu_path.split('.')[:-1]), n_production_iteration_performed, hwu_path.split('.')[-1])
                    with open(hwu_path_for_this_iteration,'w') as f:
                        f.write(tmp_obs.format_to_HwU())
                    if self.verbosity>0: logger.info("A total of %d histograms are output to file '%s' for iteration #%d. They can be rendered using the madgraph/various/histograms.py script."%(
                        len(tmp_obs),hwu_path_for_this_iteration, n_production_iteration_performed)
                    )

                    if self.verbosity>0: logger.info("Iteration #%d completed."%n_production_iteration_performed)
                    logger.info("Production results after iteration #%d/%d: %.5g +/- %.3g (n_events = %d, n_samples = %d)"%(
                        n_production_iteration_performed,n_iterations_production,
                        tmp_obs[0].histogram.bins[1].integral,
                        (tmp_obs[0].histogram.bins[1].variance**0.5),
                        tmp_obs[0].histogram.bins[1].n_entries,
                        tmp_obs[0].histogram.n_total_samples
                    ))

                # Clean up observables by removing very small weights for instance, Also apply the normalisation factor
                self.observables.clean_up(normalisation_factor=(1./n_production_iteration_performed))

                logger.info("Final production results: %.5g +/- %.3g (n_events = %d, n_samples = %d)"%(
                    self.observables[0].histogram.bins[1].integral,
                    self.observables[0].histogram.bins[1].variance**0.5,
                    self.observables[0].histogram.bins[1].n_entries,
                    self.observables[0].histogram.n_total_samples
                ))
                # Consistency check (of course will correspond to the inclusive only if range sufficiently inclusive)
                # logger.info("Final production results from ptj plot: %.5g +/- %.3g (n_events = %d, n_samples = %d)"%(
                #     self.observables[1].histogram.norm()[0],
                #     self.observables[1].histogram.norm()[1],
                #     self.observables[1].histogram.n_total_entries,
                #     self.observables[1].histogram.n_total_samples
                # ))
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
                logger.info("A total of %d histograms are output to file '%s'. They can be rendered using the madgraph/various/histograms.py script."%(
                    len(self.observables),hwu_path)
                )

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
    parser.add_argument('--min_pt', '-min_pt', dest='min_pt', type=float, default=10.,
                        help='Minimum pt imposed to the sum of final state jets (default: %(default)s).')
    parser.add_argument('--max_pt', '-max_pt', dest='max_pt', type=float, default=2000.,
                        help='Maximum pt imposed to the sum of final state jets (default: %(default)s).')
    parser.add_argument('--max_x1', '-max_x1', dest='max_x1', type=float, default=1.,
                        help='Maximum value of Bjorken x1 (default: %(default)s).')
    parser.add_argument('--max_x2', '-max_x2', dest='max_x2', type=float, default=1.,
                        help='Maximum value of Bjorken x2 (default: %(default)s).')
    parser.add_argument('--resolution_soft', '-rs', dest='resolution_soft', type=float, default=0.,
                        help='Default soft resolution parameter (default: %(default)s).')
    parser.add_argument('--resolution_coll', '-rc', dest='resolution_coll', type=float, default=0.,
                        help='Default collinear resolution parameter (default: %(default)s).')
    parser.add_argument('--separate_h_cube_errors', '-she', dest='separate_h_cube_errors', action='store_true', default=False,
                        help='Separately keep track of central values and errors for each h cube (default: disabled).')
    parser.add_argument('--no_checkoints', '-ncp', dest='use_checkpoints', action='store_false', default=True,
                        help='Disable the checkpointing system after each batch (default: enabled).')
    parser.add_argument('--n_evals_training', '-net', dest='n_evals_training', type=int, default=10000,
                        help='Number of sample points for training (default: %(default)s).')
    parser.add_argument('--n_evals_production', '-nep', dest='n_evals_production', type=int, default=10000,
                        help='Number of sample points for production (default: %(default)s).')
    parser.add_argument('--n_iterations_training', '-nit', dest='n_iterations_training', type=int, default=10,
                        help='Number of iterations for training VEGAS grids (default: %(default)s).')
    parser.add_argument('--n_iterations_production', '-nip', dest='n_iterations_production', type=int, default=3,
                        help='Number of iterations for production (default: %(default)s).')
    parser.add_argument('--n_cores', '-c', dest='n_cores', type=int, default=multiprocessing.cpu_count(),
                        help='Specify number of cores for the parallelisation (default: %(default)s).')
    parser.add_argument('--observables', '-obs', dest='observables', type=str, nargs='*', default=['ALL'],
                        help='List observables to consider. Use keyword "ALL" to select them all (default: %(default)s).')
    parser.add_argument('--tag', '-tag', dest='tag', type=str, default='tqq', choices=('st','tqg','s','tqq','u', 'all_trees'),
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
        tag = args.tag,
        resolution_coll = args.resolution_coll,
        resolution_soft = args.resolution_soft,
        separate_h_cube_errors = args.separate_h_cube_errors,
        min_pt = args.min_pt,
        max_pt = args.max_pt,
        max_x1 = args.max_x1,
        max_x2 = args.max_x2,
    )
    
    lu_dy.integrate(
        observables=args.observables,
        seed = args.seed,
        use_checkpoints = args.use_checkpoints,
        vegas_grid = args.vegas_grid,
        nhcube_batch = args.batch_size,
        hwu_path = args.hwu_path,
        n_evals_training = args.n_evals_training,
        n_iterations_training = args.n_iterations_training,
        n_iterations_production = args.n_iterations_production,
        n_evals_production = args.n_evals_production
    )
