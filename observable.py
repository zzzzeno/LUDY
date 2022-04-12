from multiprocessing.managers import BaseManager
from multiprocessing import Value

import ctypes
import vectors as vec
import logging
logger = logging.getLogger("Observable")

class ObservableError(Exception):
    pass

class Bin(object):

    def __init__(self, n_entries=0, integral=0., variance=0.):
        self.n_entries = n_entries
        self.integral = integral
        self.variance = variance

class Timings(object):
    def __init__(self, cpp=0., observables=0.):
        self.cpp = cpp
        self.observables = observables
    
    def increment(self, quantity, value):
        if quantity=='c++':
            self.cpp += value
        elif quantity=='observables':
            self.observables += value
    
    def to_dict(self):
        return {'c++':self.cpp, 'observables': self.observables}

class Histogram(object):

    def __init__(self, title, min_value, max_value, n_bins, histogram_type=None, x_axis='log', y_axis='log', **opts):
        self.title = title
        self.min = min_value
        self.max = max_value
        self.n_bins = n_bins
        self.type = histogram_type
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.n_total_events = 0
        self.n_events_rejected = 0
        self.bin_width = (self.max-self.min)/self.n_bins        
        self.bins = [Bin() for _ in range(self.n_bins)]

    def get_bins(self):
        return self.bins

    def place(self, x):
        return int((x-self.min)//self.bin_width)

    def add_weights(self, xs_and_wgts):

        wgts = {}
        for x, wgt in xs_and_wgts:
            bin_index = self.place(x)
            if bin_index in wgts:
                wgts[bin_index][0] += 1
                wgts[bin_index][1] += wgt
            else:
                wgts[bin_index] = [1,wgt]

        for bin_index, (n_new_entries, total_wgt)  in wgts.items():
            self.bins[bin_index].n_entries += n_new_entries
            self.bins[bin_index].integral += total_wgt
            self.bins[bin_index].variance += total_wgt**2

    def accumulate_histograms(self, histograms):

        for i_bin, a_bin in enumerate(self.bins):
            a_bin.integral += sum(h.bins[i_bin].integral for h in histograms.values())
            a_bin.variance += sum(
                (h.bins[i_bin].variance*h.n_total_events - h.bins[i_bin].integral**2)/(h.n_total_events-1) 
                for h in histograms.values()
            )
            a_bin.n_entries += sum(h.bins[i_bin].n_entries for h in histograms.values())

    def format_to_HwU(self):
        
        #res = ['##& xmin & xmax & central value & dy']
        res = ['<histogram> %(n_bins)d "%(title)s |X_AXIS@%(x_axis_type)s |Y_AXIS@%(y_axis_type)s |TYPE@%(histo_type)s"'%{
            'n_bins': self.n_bins, 'title': self.title, 'x_axis_type' : self.x_axis.upper(), 'y_axis_type' : self.y_axis.upper(), 'histo_type' : self.type.upper()
        },]
        for i_b, b in enumerate(self.bins):
            res.append('%.5e %.5e %.10e %.10e'%(
                self.min+i_b*self.bin_width, self.min+(i_b+1)*self.bin_width, b.integral, b.variance**0.5
            ))
        res.append("<\\histogram>")
        return '\n'.join(res)

class HistogramCollection(dict):
    def __init__(self, *args, **opts):
        pass

    def add_histogram(self, key, histo):
        self[key]= histo

class LUDYManager(BaseManager):
    pass

LUDYManager.register('Timings', Timings)
LUDYManager.register('dict', dict)
LUDYManager.register('list', list)

class Jet(object):

    def __init__(self, p, flavour):
        self.p = vec.LorentzVector(p)
        self.flavour = flavour
    
    def __str__(self):
        return "(E=%.5g, p_x=%.5g, p_y=%.5g, p_z=%.5g)@flav=%s"%(self.p[0],self.p[1],self.p[2],self.p[3],self.flavour)

class JetList(list):
    
    def __init__(self, *args, **opts):
        super(JetList, self).__init__(*args, **opts)


class Event(object):

    def __init__(self, initial_state_jets, final_state_jets, E_com, weight):
        self.initial_state_jets = initial_state_jets
        self.final_state_jets = final_state_jets
        self.weight = weight
        self.E_com = E_com
    
    def __str__(self):
        res = ["Event with %d initial-state jets and %d on final-state jet for a collision at E_com=%.4e."%(len(self.initial_state_jets),len(self.final_state_jets),self.E_com)]
        res.append("Initial state jets : %s"%(' | '.join(str(j) for j in self.initial_state_jets)))
        res.append("Final state jets   : %s"%(' | '.join(str(j) for j in self.final_state_jets)))
        res.append("Weight             : %.16e"%self.weight)
        return '\n'.join(res)

    def compute_derived_quantities(self):
        """ Compute derived quantities once here that may appear in multiple plots, for instance x_1 / x_2."""
        #TODO
        pass

# Keep events from a single sample point grouped so as to properly account for correlations
class EventGroup(list):
    def __init__(self, *args, h_cube=None, **opts):
        self.h_cube =h_cube
        super(EventGroup, self).__init__(*args, **opts)

    def compute_derived_quantities(self, *args, **opts):
        for e in self:
            e.compute_derived_quantities(*args, **opts)

    def __str__(self):
        res = ["Event group for h_cube #%d containing %d events:"%(self.h_cube, len(self))]
        for i_e, e in enumerate(self):
            res.extend(['     %s'%l if i_l > 0 else '#%-3d %s'%(i_e,l) for i_l, l in enumerate(str(e).split('\n'))])
        return '\n'.join(res)

class Observable(object):

    def __init__(self, *args, **opts):
        self.histogram_args = args
        self.histogram_opts = opts
        self.histogram = Histogram(*self.histogram_args, **self.histogram_opts)
        
        self.histos_current_iteration = {}
        
        super(Observable, self).__init__()

    def accumulate_event_group(self, event_group):
        if self.histos_current_iteration is None:
            raise ObservableError("You must call 'start_iteration()' before starting to accumulate events.")
        if event_group.h_cube not in self.histos_current_iteration.keys():
            self.histos_current_iteration[event_group.h_cube] = Histogram(*self.histogram_args, **self.histogram_opts)
        histogram = self.histos_current_iteration.get(event_group.h_cube)
        histogram.add_weights( [(self(e), e.weight) for e in event_group] )

    def __call__(self, event):
        raise NotImplementedError("This function must be implemented by the daughter class.")

    def finalize_iteration(self, worker_observable=None):

        if worker_observable is not None:
            histos_to_merge = worker_observable.histos_current_iteration
        else:
            histos_to_merge = self.histos_current_iteration

        if len(histos_to_merge)==0:
            return
        
        self.histogram.accumulate_histograms(histos_to_merge)
        histos_to_merge.clear()

    def format_to_HwU(self):
        return self.histogram.format_to_HwU()

class ObservableList(list):
    
    def __init__(self, *args, **opts):
        super(ObservableList, self).__init__(*args, **opts)
    
    def start_iteration(self, *args, **opts):
        for o in self:
            o.start_iteration(*args, **opts)

    def accumulate_event_group(self, *args, **opts):
        for o in self:
            o.accumulate_event_group(*args, **opts)

    def finalize_iteration(self, *args, worker_observable=None, **opts):
        if worker_observable is None:
            for o in self:
                o.finalize_iteration(*args, **opts)
        else:
            for o, w_o in zip(self, worker_observable):
                o.finalize_iteration(*args, worker_observable = w_o, **opts)

    def format_to_HwU(self):
        res = ["##& xmin & xmax & central value & dy",'']
        for o in self:
            res.append(o.format_to_HwU())
        return '\n'.join(res)

class CrossSection(Observable):

    def __init__(self, title='CrossSection', histogram_type='NLO', **opts):
        super(CrossSection, self).__init__(title, min_value=0., max_value=3., n_bins=3, histogram_type=histogram_type, x_axis='lin', y_axis='lin', **opts)

    def __call__(self, event):
        return 1.5

class x1(Observable):

    def __init__(self, title='x1', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='log', **opts):
        super(x1, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        # TODO
        return 1.0
