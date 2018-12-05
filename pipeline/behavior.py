'''
Schema of session information.
'''
import datajoint as dj
from pipeline import acquisition

schema = dj.schema('ttngu207_behavior',locals())

        
@schema
class TrialType(dj.Lookup):
    definition = """
    type_name: varchar(64)
    ---
    
    """
    contents = [
            ['Hit'],
            ['Miss'],
            ['CR'],
            ['FA'],
            ['Stimtrials']
            ]

@schema
class TrialSet(dj.Manual):
    definition = """
    -> acquisition.Session
    ---
    number_of_trials: int
    time_unit: enum('millisecond','second','minute','hour','day')
    """
    class Trial(dj.Part):
        definition = """
        -> master
        trial_idx: int
        ---
        -> TrialType
        pole_trial_condition: enum('Go','NoGo')  # string indicating whether the pole was presented in a ‘Go’ or ‘Nogo’ location
        pole_position: float                     # is the location of the pole along the anteroposterior axis of the animal in microsteps (0.0992 microns / microstep)
        pole_in_time: float                      # is the start of sample period for each trial (e.g. the onset of pole motion towards the exploration area), in units of seconds, relative to trialStartTimes
        pole_out_time: float                     # is the end of the sample period (e.g. the onset of pole motion away from the exploration area)
        lick_time: longblob                      # is an array of times of when the mouse’s tongue initiates contact with the spout
        
        start_sample: int       # the index of the starting sample of this trial, with respect to the starting of this session (at 0th sample point) - this way, time will be derived from the sampling rate of a particular recording downstream
        end_sample: int         # the index of the ending sample of this trial, with respect to the starting of this session (at 0th sample point) - this way, time will be derived from the sampling rate of a particular recording downstream
        
        """





















    


   
    
    