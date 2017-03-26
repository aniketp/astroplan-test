from astroplan import Observer
from astroplan import FixedTarget
from astropy.time import Time
from astroplan.constraints import AtNightConstraint, AirmassConstraint, AltitudeConstraint
from astroplan import ObservingBlock
from astroplan.constraints import TimeConstraint
from astropy import units as u
from astroplan.scheduling import Transitioner
from astroplan.scheduling import Schedule
from abc import ABCMeta, abstractmethod
import copy
import numpy as np

class Scheduler(object):
    __metaclass__ = ABCMeta

    @u.quantity_input(gap_time=u.second, time_resolution=u.second)
    def __init__(self, constraints, observer, transitioner=None,
                 gap_time=5*u.min, time_resolution=20*u.second):
       
        self.constraints = constraints
        self.observer = observer
        self.transitioner = transitioner
        if not isinstance(self.transitioner, Transitioner):
            raise ValueError("A Transitioner is required")
        self.gap_time = gap_time
        self.time_resolution = time_resolution

    def __call__(self, blocks, schedule):
        
        self.schedule = schedule
        self.schedule.observer = self.observer
        
        copied_blocks = [copy.copy(block) for block in blocks]
        schedule = self._make_schedule(copied_blocks)
        return schedule

    @abstractmethod
    def _make_schedule(self, blocks):
        
        raise NotImplementedError
        return schedule

    @classmethod
    @u.quantity_input(duration=u.second)
    def from_timespan(cls, center_time, duration, **kwargs):
        start_time = center_time - duration / 2.
        end_time = center_time + duration / 2.
        return cls(start_time, end_time, **kwargs)

class SequentialScheduler(Scheduler):
    def __init__(self, *args, **kwargs):
        super(SequentialScheduler, self).__init__(*args, **kwargs)
    
    @profile
    def _make_schedule(self, blocks):
        pre_filled = np.array([[block.start_time, block.end_time] for
                      block in self.schedule.scheduled_blocks])
        if len(pre_filled) == 0:
            a = self.schedule.start_time
            filled_times = Time([a - 1*u.hour, a - 1*u.hour,
                               a - 1*u.minute, a - 1*u.minute])
            pre_filled = filled_times.reshape((2, 2))
        else:
            filled_times = Time(pre_filled.flatten())
            pre_filled = filled_times.reshape((int(len(filled_times)/2), 2))
        for b in blocks:
            if b.constraints is None:
                b._all_constraints = self.constraints
            else:
                b._all_constraints = self.constraints + b.constraints
            
            if b._all_constraints is None:
                b._all_constraints = [AltitudeConstraint(min=0 * u.deg)]
                b.constraints = [AltitudeConstraint(min=0 * u.deg)]
            elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
                b._all_constraints.append(AltitudeConstraint(min=0 * u.deg))
                if b.constraints is None:
                    b.constraints = [AltitudeConstraint(min=0 * u.deg)]
                else:
                    b.constraints.append(AltitudeConstraint(min=0 * u.deg))
            b._duration_offsets = u.Quantity([0*u.second, b.duration/2,
                                              b.duration])
            b.observer = self.observer
        current_time = self.schedule.start_time
        while (len(blocks) > 0) and (current_time < self.schedule.end_time):
           
            block_transitions = []
            block_constraint_results = []
            for b in blocks:
                
                if len(self.schedule.observing_blocks) > 0:
                    trans = self.transitioner(self.schedule.observing_blocks[-1], b, current_time,
                                              self.observer)
                else:
                    trans = None
                block_transitions.append(trans)
                transition_time = 0*u.second if trans is None else trans.duration

                times = current_time + transition_time + b._duration_offsets

                
                if (any((current_time < filled_times) & (filled_times < times[2])) or
                        any(abs(pre_filled.T[0]-current_time) < 1*u.second)):
                    block_constraint_results.append(0)

                else:
                    constraint_res = []
                    for constraint in b._all_constraints:
                        constraint_res.append(constraint(self.observer, [b.target],
                                                         times))
                    
                    block_constraint_results.append(np.prod(constraint_res))

            
            bestblock_idx = np.argmax(block_constraint_results)

            if block_constraint_results[bestblock_idx] == 0.:
                
                current_time += self.gap_time
            else:
                
                trans = block_transitions.pop(bestblock_idx)
                if trans is not None:
                    self.schedule.insert_slot(trans.start_time, trans)
                    current_time += trans.duration

                
                newb = blocks.pop(bestblock_idx)
                newb.start_time = current_time
                current_time += newb.duration
                newb.end_time = current_time
                newb.constraints_value = block_constraint_results[bestblock_idx]

                self.schedule.insert_slot(newb.start_time, newb)

        return self.schedule

apo = Observer.at_site('apo')

deneb = FixedTarget.from_name('Deneb')
m13 = FixedTarget.from_name('M13')

noon_before = Time('2016-07-06 19:00')
noon_after = Time('2016-07-07 19:00')

global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False), AtNightConstraint.twilight_civil()]

read_out = 20 * u.second

deneb_exp = 60*u.second
m13_exp = 100*u.second
n = 16
blocks = []

half_night_start = Time('2016-07-07 02:00')
half_night_end = Time('2016-07-07 08:00')
first_half_night = TimeConstraint(half_night_start, half_night_end)

for priority, bandpass in enumerate(['B', 'G', 'R']):
    
    b = ObservingBlock.from_exposures(deneb, priority, deneb_exp, n, read_out, configuration = {'filter': bandpass}, constraints = [first_half_night])
    blocks.append(b)

    b = ObservingBlock.from_exposures(m13, priority, m13_exp, n, read_out, configuration = {'filter': bandpass}, constraints = [first_half_night])
    blocks.append(b)

slew_rate = .8*u.deg/u.second
transitioner = Transitioner(slew_rate, {'filter':{('B','G'): 10*u.second, ('G','R'): 10*u.second, 'default': 30*u.second}})

seq_scheduler = SequentialScheduler(constraints = global_constraints, observer = apo, transitioner = transitioner)

sequential_schedule = Schedule(noon_before, noon_after)

print(seq_scheduler(blocks, sequential_schedule))
print(sequential_schedule.to_table())
