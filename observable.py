import ctypes
import vectors as vec
import logging
import copy
import math
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
    
    def aggregate_timings(self, other):
        self.cpp += other.cpp
        self.observables += other.observables

    def increment(self, quantity, value):
        if quantity=='c++':
            self.cpp += value
        elif quantity=='observables':
            self.observables += value
    
    def to_dict(self):
        return {'c++':self.cpp, 'observables': self.observables}

class Histogram(object):

    def __init__(self, title, min_value, max_value, n_bins, histogram_type=None, x_axis='lin', y_axis='log', **opts):
        self.title = title
        self.min = min_value
        self.max = max_value
        self.n_bins = n_bins
        self.type = histogram_type
        self.x_axis = x_axis
        if self.x_axis != 'lin':
            raise ObservableError("For now HwU histograms only support linear x-axis. Simply fill in log10 of the observable on the x-axis to emulate what you want.")
        self.y_axis = y_axis
        self.n_total_samples = 0
        self.n_total_entries = 0
        self.bin_width = (self.max-self.min)/self.n_bins        
        self.bins = [Bin() for _ in range(self.n_bins)]

    def norm(self):
        return sum(b.integral for b in self.bins), sum(b.variance for b in self.bins)**0.5

    def get_bins(self):
        return self.bins

    def place(self, x):
        if x >= self.max or x < self.min:
            return None
        return int((x-self.min)//self.bin_width)

    def clean_up(self, *args, normalisation_factor=1., threshold=1.e-15, **opts):
        abs_inclusive = sum(abs(b.integral) for b in self.bins)
        if abs_inclusive == 0.:
            return
        for b in self.bins:
            if abs(b.integral/abs_inclusive)<threshold:
                b.integral = 0.
                b.variance = 0.
            
            b.integral *= normalisation_factor
            b.variance *= normalisation_factor**2

    def add_weights(self, xs_and_wgts):

        wgts = {}
        self.n_total_samples += 1
        self.n_total_entries += len(xs_and_wgts)
        for x, wgt in xs_and_wgts:
            bin_index = self.place(x)
            if bin_index is None:
                continue
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
        
        self.n_total_samples += sum(h.n_total_samples for h in histograms.values())
        self.n_total_entries += sum(h.n_total_entries for h in histograms.values())
        for i_bin, a_bin in enumerate(self.bins):
            a_bin.integral += sum(h.bins[i_bin].integral for h in histograms.values())
            a_bin.variance += sum(
                (h.bins[i_bin].variance*h.n_total_samples - h.bins[i_bin].integral**2)/(h.n_total_samples-1) 
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
                self.min+i_b*self.bin_width, self.min+(i_b+1)*self.bin_width, b.integral, (b.variance*(self.n_total_entries/b.n_entries if b.n_entries > 0 else 1.))**0.5
            ))
        res.append("<\\histogram>")
        return '\n'.join(res)

class HistogramCollection(dict):
    def __init__(self, *args, **opts):
        pass

    def add_histogram(self, key, histo):
        self[key]= histo

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
        self.derived_quantities_computed = False
    
    def __str__(self):
        res = ["Event with %d initial-state jets and %d on final-state jet for a collision at E_com=%.4e."%(len(self.initial_state_jets),len(self.final_state_jets),self.E_com)]
        res.append("Initial state jets : %s"%(' | '.join(str(j) for j in self.initial_state_jets)))
        res.append("Final state jets   : %s"%(' | '.join(str(j) for j in self.final_state_jets)))
        res.append("Weight             : %.16e"%self.weight)
        if self.derived_quantities_computed:
            res.append(self.str_derived_quantities())
        return '\n'.join(res)

    def str_derived_quantities(self):
        res = ['Derived quantities :']
        res.append('  > x1             : %.16e'%self.x1)
        res.append('  > x2             : %.16e'%self.x2)
        return '\n'.join(res)

    def compute_derived_quantities(self):
        """ Compute derived quantities once here that may appear in multiple plots, for instance x_1 / x_2."""
        self.derived_quantities_computed = True
        # Store x1 and x2 for later use
        sum_isr_jets_p = sum(j.p for j in self.initial_state_jets)
        x1px2 = 2*(sum_isr_jets_p[0])/self.E_com
        x1mx2 = 2*(sum_isr_jets_p[3])/self.E_com
        self.x1=(x1px2+x1mx2)/2
        self.x2=(x1px2-x1mx2)/2

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

    def reset(self,*args,**opts):
        self.histogram = Histogram(*self.histogram_args, **self.histogram_opts)
        self.histos_current_iteration = {}

    def accumulate_event_group(self, event_group):
        if self.histos_current_iteration is None:
            raise ObservableError("You must call 'start_iteration()' before starting to accumulate events.")
        if event_group.h_cube not in self.histos_current_iteration.keys():
            self.histos_current_iteration[event_group.h_cube] = Histogram(*self.histogram_args, **self.histogram_opts)
        histogram = self.histos_current_iteration.get(event_group.h_cube)
        histogram.add_weights( sum([ self(e) for e in event_group ],[]) )

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

    def clean_up(self, *args, **opts):
        self.histogram.clean_up(*args, **opts)

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

    def reset(self, *args, **opts):
        for o in self:
            o.reset(*args, **opts)

    def finalize_iteration(self, *args, worker_observable=None, **opts):
        if worker_observable is None:
            for o in self:
                o.finalize_iteration(*args, **opts)
        else:
            for o, w_o in zip(self, worker_observable):
                o.finalize_iteration(*args, worker_observable = w_o, **opts)

    def clean_up(self, *args, **opts):
        """ Clean up histograms, for instance by removing small weights."""
        for o in self:
            o.clean_up(*args, **opts)

    @classmethod
    def merge(cls, list_of_obs_list):
        
        res=copy.deepcopy(list_of_obs_list[0])
        for ol in list_of_obs_list[1:]:
            for i_obs, obs in enumerate(ol):
                for k, h in obs.histos_current_iteration.items():
                    if k not in res[i_obs].histos_current_iteration:
                        res[i_obs].histos_current_iteration[k] = copy.deepcopy(h)
                    else:
                        res[i_obs].histos_current_iteration[k].n_total_samples += h.n_total_samples
                        res[i_obs].histos_current_iteration[k].n_total_entries += h.n_total_entries
                        for i_b, b in enumerate(res[i_obs].histos_current_iteration[k].bins):
                            b.n_entries += h.bins[i_b].n_entries
                            b.integral += h.bins[i_b].integral
                            b.variance += h.bins[i_b].variance

        return res

    def format_to_HwU(self):
        res = ["##& xmin & xmax & central value & dy",'']
        for o in self:
            res.append(o.format_to_HwU())
        return '\n'.join(res)

class CrossSection(Observable):

    def __init__(self, title='CrossSection', histogram_type='NLO', **opts):
        super(CrossSection, self).__init__(title, min_value=0., max_value=3., n_bins=3, histogram_type=histogram_type, x_axis='lin', y_axis='lin', **opts)

    def __call__(self, event):
        return [(1.5,event.weight),]

class ptj(Observable):

    def __init__(self, title='ptj', min_value=0., max_value=5000., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='log', **opts):
        super(ptj, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        # The first and unique final state jets should always be the leading jet; we could verify this if need be
        return [(event.final_state_jets[0].p.pt(), event.weight),] if len(event.final_state_jets)>0 else []

class log10ptj(Observable):

    def __init__(self, title='log10ptj', min_value=-1., max_value=4., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        super(log10ptj, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        return [( math.log10(event.final_state_jets[0].p.pt()), event.weight ),] if len(event.final_state_jets)>0 else []

class x1(Observable):

    def __init__(self, title='x1', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        super(x1, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        return [(event.x1, event.weight),]

class x2(Observable):

    def __init__(self, title='x2', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        super(x2, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        return [(event.x2, event.weight),]

class z(Observable):

    def __init__(self, title='z', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        super(z, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        return [(event.x1*event.x2, event.weight),]

class log10z(Observable):

    def __init__(self, title='z', min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        super(log10z, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        return [( math.log10(event.x1*event.x2), event.weight),]

class x1FixedFlav(Observable):

    def __init__(self, flavour=0, take_abs=True, min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        title = 'x1%sFlavour%s'%('Abs' if take_abs else '','m%d'%abs(flavour) if flavour<0 else '0' if flavour==0 else 'p%d'%flavour)
        self.take_abs = take_abs
        self.flavour = flavour
        super(x1FixedFlav, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        #TODO check if event.initial_state_jets[0].flavour really always correspond to event.x1 indeed!
        return [(event.x1, event.weight),] if (abs(event.initial_state_jets[0].flavour) if self.take_abs else event.initial_state_jets[0].flavour)==self.flavour else []

class x2FixedFlav(Observable):

    def __init__(self, flavour=0, take_abs=True, min_value=0., max_value=1., n_bins=100, histogram_type='NLO', x_axis='lin', y_axis='lin', **opts):
        title = 'x2%sFlavour%s'%('Abs' if take_abs else '','m%d'%abs(flavour) if flavour<0 else '0' if flavour==0 else 'p%d'%flavour)
        self.take_abs = take_abs
        self.flavour = flavour
        super(x2FixedFlav, self).__init__(title, min_value, max_value, n_bins, histogram_type=histogram_type, x_axis=x_axis, y_axis=y_axis, **opts)

    def __call__(self, event):
        #TODO check if event.initial_state_jets[1].flavour really always correspond to event.x2 indeed!
        return [(event.x2, event.weight),] if (abs(event.initial_state_jets[1].flavour) if self.take_abs else event.initial_state_jets[0].flavour)==self.flavour else []