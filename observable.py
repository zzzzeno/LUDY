from multiprocessing.managers import BaseManager
from multiprocessing import Value

import ctypes
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

    def __init__(self, p, spin):
        self.p = p
        self.spin = spin

class JetList(list):
    
    def __init__(self, *args, **opts):
        super(JetList, self).__init__(*args, **opts)


class Event(object):

    def __init__(self, initial_state_jets, final_state_jets, weight):
        self.initial_state_jets = initial_state_jets
        self.final_state_jets = final_state_jets
        self.weight = weight

# Keep events from a single sample point grouped so as to properly account for correlations
class EventGroup(list):
    def __init__(self, *args, h_cube=None, **opts):
        self.h_cube =h_cube
        super(EventGroup, self).__init__(*args, **opts)

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
        self.accumulate(histogram, event_group)

    def accumulate(self, histogram, event):
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

class CrossSection(Observable):

    def __init__(self, title='CrossSection', histogram_type='NLO', **opts):
        super(CrossSection, self).__init__(title, min_value=0., max_value=3., n_bins=3, histogram_type=histogram_type, x_axis='lin', y_axis='lin', **opts)

    def accumulate(self, histogram, event_group):
        histogram.add_weights( [(1.5, e.weight) for e in event_group] )

class x1(Observable):

    def __init__(self, beam_energy, title='x1', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='log', **opts):
        self.beam_energy = beam_energy
        super(x1, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def accumulate(self, histogram, event):
        #TODO
        pass