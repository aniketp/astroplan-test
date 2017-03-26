
Schedule containing 6 observing blocks between 2016-07-06 19:00:00.000 and 2016-07-07 19:00:00.000
     target         start time (UTC)         end time (UTC)     duration (minutes)        ra            dec         configuration  
--------------- ----------------------- ----------------------- ------------------ --------------- -------------- -----------------
            M13 2016-07-07 02:45:00.000 2016-07-07 03:17:00.000               32.0 250d25m24.4499s 36d27m40.7002s   {'filter': 'B'}
TransitionBlock 2016-07-07 03:17:00.000 2016-07-07 03:17:30.000                0.5                                ['filter:B to R']
            M13 2016-07-07 03:17:30.000 2016-07-07 03:49:30.000               32.0 250d25m24.4499s 36d27m40.7002s   {'filter': 'R'}
TransitionBlock 2016-07-07 03:49:30.000 2016-07-07 03:50:00.000                0.5                                ['filter:R to G']
            M13 2016-07-07 03:50:00.000 2016-07-07 04:22:00.000               32.0 250d25m24.4499s 36d27m40.7002s   {'filter': 'G'}
TransitionBlock 2016-07-07 04:22:00.000 2016-07-07 04:23:26.384      1.43973064423                                ['filter:G to B']
          Deneb 2016-07-07 04:23:26.384 2016-07-07 04:44:46.384      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'B'}
TransitionBlock 2016-07-07 04:44:46.384 2016-07-07 04:45:16.384                0.5                                ['filter:B to R']
          Deneb 2016-07-07 04:45:16.384 2016-07-07 05:06:36.384      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'R'}
TransitionBlock 2016-07-07 05:06:36.384 2016-07-07 05:07:06.384                0.5                                ['filter:R to G']
          Deneb 2016-07-07 05:07:06.384 2016-07-07 05:28:26.384      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'G'}
Wrote profile results to seq_test.py.lprof
Timer unit: 1e-06 s

Total time: 37.6822 s
File: seq_test.py
Function: _make_schedule at line 55

#### Sequential Scheduling

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    55                                               @profile
    56                                               def _make_schedule(self, blocks):
    57         1            5      5.0      0.0          pre_filled = np.array([[block.start_time, block.end_time] for
    58         1           22     22.0      0.0                        block in self.schedule.scheduled_blocks])
    59         1            4      4.0      0.0          if len(pre_filled) == 0:
    60         1            2      2.0      0.0              a = self.schedule.start_time
    61         1         3543   3543.0      0.0              filled_times = Time([a - 1*u.hour, a - 1*u.hour,
    62         1         3369   3369.0      0.0                                 a - 1*u.minute, a - 1*u.minute])
    63         1          141    141.0      0.0              pre_filled = filled_times.reshape((2, 2))
    64                                                   else:
    65                                                       filled_times = Time(pre_filled.flatten())
    66                                                       pre_filled = filled_times.reshape((int(len(filled_times)/2), 2))
    67         7           14      2.0      0.0          for b in blocks:
    68         6           12      2.0      0.0              if b.constraints is None:
    69                                                           b._all_constraints = self.constraints
    70                                                       else:
    71         6           16      2.7      0.0                  b._all_constraints = self.constraints + b.constraints
    72                                                       
    73         6            8      1.3      0.0              if b._all_constraints is None:
    74                                                           b._all_constraints = [AltitudeConstraint(min=0 * u.deg)]
    75                                                           b.constraints = [AltitudeConstraint(min=0 * u.deg)]
    76         6           40      6.7      0.0              elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
    77                                                           b._all_constraints.append(AltitudeConstraint(min=0 * u.deg))
    78                                                           if b.constraints is None:
    79                                                               b.constraints = [AltitudeConstraint(min=0 * u.deg)]
    80                                                           else:
    81                                                               b.constraints.append(AltitudeConstraint(min=0 * u.deg))
    82         6         1000    166.7      0.0              b._duration_offsets = u.Quantity([0*u.second, b.duration/2,
    83         6         1195    199.2      0.0                                                b.duration])
    84         6           17      2.8      0.0              b.observer = self.observer
    85         1            2      2.0      0.0          current_time = self.schedule.start_time
    86       100         3877     38.8      0.0          while (len(blocks) > 0) and (current_time < self.schedule.end_time):
    87                                                      
    88        99          180      1.8      0.0              block_transitions = []
    89        99          150      1.5      0.0              block_constraint_results = []
    90       678         1145      1.7      0.0              for b in blocks:
    91                                                           
    92       579         3895      6.7      0.0                  if len(self.schedule.observing_blocks) > 0:
    93        15      1171904  78126.9      3.1                      trans = self.transitioner(self.schedule.observing_blocks[-1], b, current_time, self.observer)
    94                                                           else:
    95       564          637      1.1      0.0                      trans = None
    96       579          843      1.5      0.0                  block_transitions.append(trans)
    97       579        28459     49.2      0.1                  transition_time = 0*u.second if trans is None else trans.duration
    98                                           
    99       579       876907   1514.5      2.3                  times = current_time + transition_time + b._duration_offsets
   100                                           
   101                                                           
   102       579       107218    185.2      0.3                  if (any((current_time < filled_times) & (filled_times < times[2])) or
   103       579       674889   1165.6      1.8                          any(abs(pre_filled.T[0]-current_time) < 1*u.second)):
   104                                                               block_constraint_results.append(0)
   105                                           
   106                                                           else:
   107       579         1924      3.3      0.0                      constraint_res = []
   108      2316         4917      2.1      0.0                      for constraint in b._all_constraints:
   109      1737     34562650  19897.9     91.7                          constraint_res.append(constraint(self.observer, [b.target], times))
   110                                                               
   111       579        16248     28.1      0.0                      block_constraint_results.append(np.prod(constraint_res))
   112                                           
   113                                                       
   114        99         2160     21.8      0.0              bestblock_idx = np.argmax(block_constraint_results)
   115                                           
   116        99          201      2.0      0.0              if block_constraint_results[bestblock_idx] == 0.:
   117                                                           
   118        93        70354    756.5      0.2                  current_time += self.gap_time
   119                                                       else:
   120                                                           
   121         6           14      2.3      0.0                  trans = block_transitions.pop(bestblock_idx)
   122         6            6      1.0      0.0                  if trans is not None:
   123         5        61009  12201.8      0.2                      self.schedule.insert_slot(trans.start_time, trans)
   124         5         3843    768.6      0.0                      current_time += trans.duration
   125                                           
   126                                                           
   127         6           20      3.3      0.0                  newb = blocks.pop(bestblock_idx)
   128         6           12      2.0      0.0                  newb.start_time = current_time
   129         6         4295    715.8      0.0                  current_time += newb.duration
   130         6           12      2.0      0.0                  newb.end_time = current_time
   131         6           10      1.7      0.0                  newb.constraints_value = block_constraint_results[bestblock_idx]
   132                                           
   133         6        74989  12498.2      0.2                  self.schedule.insert_slot(newb.start_time, newb)
   134                                           
   135         1            1      1.0      0.0          return self.schedule

#### Scheduler call

Total time: 45.4981 s
File: seq_test.py
Function: __call__ at line 29

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    29                                               @profile
    30                                               def __call__(self, blocks, schedule):
    31                                                   
    32         1            8      8.0      0.0          self.schedule = schedule
    33         1            4      4.0      0.0          self.schedule.observer = self.observer
    34                                                   
    35         1          537    537.0      0.0          copied_blocks = [copy.copy(block) for block in blocks]
    36         1     45497577 45497577.0    100.0          schedule = self._make_schedule(copied_blocks)
    37         1            2      2.0      0.0          return schedule
    
 
#### Alt-Az

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   146                                               @profile
   147                                               def altaz(self, time, target=None, obswl=None, grid=True):
   148                                                  
   149         1            2      2.0      0.1          if target is not None:
   150                                                       time, target = self._preprocess_inputs(time, target, grid)
   151                                           
   152         1         3565   3565.0     99.9          altaz_frame = AltAz(location=self.location, obstime=time,pressure=self.pressure, obswl=obswl, temperature=self.temperature, relative_humidity=self.relative_humidity)
   153         1            1      1.0      0.0          if target is None:
   154                                                       # Return just the frame
   155         1            0      0.0      0.0              return altaz_frame
   156                                                   else:
   157                                                       return target.transform_to(altaz_frame)

  
