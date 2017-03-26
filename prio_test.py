from astroplan import Observer
from astroplan import FixedTarget
from astropy.time import Time
from astroplan.constraints import AtNightConstraint, AirmassConstraint, AltitudeConstraint
from astroplan import ObservingBlock
from astroplan.constraints import TimeConstraint
from astropy import units as u
from astroplan.scheduling import Transitioner, TransitionBlock
from astroplan.scheduling import Schedule
from abc import ABCMeta, abstractmethod
import copy
import numpy as np
from astroplan import Scorer

def stride_array(arr, window_width):
    as_strided = np.lib.stride_tricks.as_strided

    new_shape = (len(arr) - window_width + 1, window_width)

    strided_arr = as_strided(arr, new_shape, (arr.strides[0], arr.strides[0]))

    return strided_arr


def time_grid_from_range(time_range, time_resolution=0.5*u.hour):
    
    try:
        start_time, end_time = time_range
    except ValueError:
        raise ValueError("time_range should have a length of 2: lower and "
                         "upper bounds on the time sequence.")
    return Time(np.arange(start_time.jd, end_time.jd,
                          time_resolution.to(u.day).value), format='jd')

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

class PriorityScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        """

        """
        super(PriorityScheduler, self).__init__(*args, **kwargs)

    @profile
    def _make_schedule(self, blocks):

        _all_times = []
        _block_priorities = np.zeros(len(blocks))
        for i, b in enumerate(blocks):
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
            b._duration_offsets = u.Quantity([0 * u.second, b.duration / 2, b.duration])
            _block_priorities[i] = b.priority
            _all_times.append(b.duration)
            b.observer = self.observer

       

        time_resolution = self.time_resolution
        times = time_grid_from_range([self.schedule.start_time, self.schedule.end_time],
                                     time_resolution=time_resolution)
        is_open_time = np.ones(len(times), bool)
        # close times that are already filled
        pre_filled = np.array([[block.start_time, block.end_time] for
                               block in self.schedule.scheduled_blocks])
        for start_end in pre_filled:
            filled = np.where((start_end[0] < times) & (times < start_end[1]))
            is_open_time[filled[0]] = False
            is_open_time[min(filled[0]) - 1] = False
        # generate the score arrays for all of the blocks
        scorer = Scorer(blocks, self.observer, self.schedule, global_constraints=self.constraints)
        score_array = scorer.create_score_array(time_resolution)

        # Sort the list of blocks by priority
        sorted_indices = np.argsort(_block_priorities)

        unscheduled_blocks = []
        # Compute the optimal observation time in priority order
        for i in sorted_indices:
            b = blocks[i]
            # Compute possible observing times by combining object constraints
            # with the master open times mask
            constraint_scores = score_array[i]

            # Add up the applied constraints to prioritize the best blocks
            # And then remove any times that are already scheduled
            constraint_scores[is_open_time == False] = 0
            # Select the most optimal time

            # need to leave time around the Block for transitions
            if self.transitioner.instrument_reconfig_times:
                max_config_time = sum([max(value.values()) for value in
                                       self.transitioner.instrument_reconfig_times.values()])
            else:
                max_config_time = 0*u.second
            if self.transitioner.slew_rate:
                buffer_time = (160*u.deg/self.transitioner.slew_rate + max_config_time)
            else:
                buffer_time = max_config_time
            # TODO: make it so that this isn't required to prevent errors in slot creation
            total_duration = b.duration + buffer_time
            # calculate the number of time slots needed for this exposure
            _stride_by = np.int(np.ceil(float(total_duration / time_resolution)))

            # Stride the score arrays by that number
            _strided_scores = stride_array(constraint_scores, _stride_by)

           
            good = np.all(_strided_scores > 1e-5, axis=1)
            sum_scores = np.zeros(len(_strided_scores))
            sum_scores[good] = np.sum(_strided_scores[good], axis=1)

            if np.all(constraint_scores == 0) or np.all(good == False):
                # No further calculation if no times meet the constraints
                _is_scheduled = False

            else:
                # If an optimal block is available, _is_scheduled=True
                best_time_idx = np.argmax(sum_scores)
                start_time_idx = best_time_idx
                new_start_time = times[best_time_idx]
                _is_scheduled = True

            if _is_scheduled:
                # set duration such that the Block will fit in the strided array
                duration_indices = np.int(np.ceil(float(b.duration / time_resolution)))
                b.duration = duration_indices * time_resolution
                # add 1 second to the start time to allow for scheduling at the start of a slot
                slot_index = [q for q, slot in enumerate(self.schedule.slots)
                              if slot.start < new_start_time + 1*u.second < slot.end][0]
                slots_before = self.schedule.slots[:slot_index]
                slots_after = self.schedule.slots[slot_index + 1:]
                # this has to remake transitions between already existing ObservingBlocks
                if slots_before:
                    if isinstance(self.schedule.slots[slot_index - 1].block, ObservingBlock):
                        # make a transition object after the previous ObservingBlock
                        tb = self.transitioner(self.schedule.slots[slot_index - 1].block, b,
                                               self.schedule.slots[slot_index - 1].end, self.observer)
                        times_indices = np.int(np.ceil(float(tb.duration / time_resolution)))
                        tb.duration = times_indices * time_resolution
                        start_idx = self.schedule.slots[slot_index - 1].block.end_idx
                        end_idx = times_indices + start_idx
                        # this may make some OBs get sub-optimal scheduling, but it closes gaps
                        # TODO: determine a reasonable range inside which it gets shifted
                        if (new_start_time - tb.start_time < tb.duration or
                                abs(new_start_time - tb.end_time) < self.gap_time):
                            new_start_time = tb.end_time
                            start_time_idx = end_idx
                        self.schedule.insert_slot(tb.start_time, tb)
                        is_open_time[start_idx: end_idx] = False
                        slot_index += 1
                        # Remove times from the master time list (copied in later code blocks)
                    elif isinstance(self.schedule.slots[slot_index - 1].block, TransitionBlock):
                        # change the existing TransitionBlock to what it needs to be now
                        tb = self.transitioner(self.schedule.slots[slot_index - 2].block, b,
                                               self.schedule.slots[slot_index - 2].end, self.observer)
                        times_indices = np.int(np.ceil(float(tb.duration / time_resolution)))
                        tb.duration = times_indices * time_resolution
                        start_idx = self.schedule.slots[slot_index - 2].block.end_idx
                        end_idx = times_indices + start_idx
                        self.schedule.change_slot_block(slot_index - 1, new_block=tb)
                        if (new_start_time - tb.start_time < tb.duration or
                                abs(new_start_time - tb.end_time) < self.gap_time):
                            new_start_time = tb.end_time
                            start_time_idx = end_idx
                        is_open_time[start_idx: end_idx] = False
                end_time_idx = duration_indices + start_time_idx

                if slots_after:
                    if isinstance(self.schedule.slots[slot_index + 1].block, ObservingBlock):
                        # make a transition object after the new ObservingBlock
                        tb = self.transitioner(b, self.schedule.slots[slot_index + 1].block,
                                               new_start_time + b.duration, self.observer)
                        times_indices = np.int(np.ceil(float(tb.duration / time_resolution)))
                        tb.duration = times_indices * time_resolution
                        self.schedule.insert_slot(tb.start_time, tb)
                        start_idx = end_time_idx
                        end_idx = start_idx + times_indices
                        is_open_time[start_idx: end_idx] = False

                # now assign the block itself times and add it to the schedule
                b.constraints = b._all_constraints
                b.end_idx = end_time_idx
                self.schedule.insert_slot(new_start_time, b)
                is_open_time[start_time_idx: end_time_idx] = False

            else:
                print("could not schedule", b.target.name)
                unscheduled_blocks.append(b)
                continue

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

pri_scheduler = PriorityScheduler(constraints = global_constraints, observer = apo, transitioner = transitioner)

priority_schedule = Schedule(noon_before, noon_after)

print(pri_scheduler(blocks, priority_schedule))
print(priority_schedule.to_table())
