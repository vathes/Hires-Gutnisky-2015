'''
Schema of aquisition information.
'''
import datajoint as dj
from pipeline import reference, subject

schema = dj.schema('ttngu207_acquisition',locals())

@schema
class ExperimentType(dj.Lookup):
    definition = """
    experiment_type_name: varchar(64)
    """
    contents = [
        ['behavior'], ['extracelluar'], ['photostim']
    ]

@schema
class PhotoStim(dj.Manual):
    definition = """
    photo_stim_id: int
    ---
    photo_stim_wavelength: int
    photo_stim_method: enum('fiber', 'laser')
    -> reference.BrainLocation.proj(photo_stim_location="brain_location")
    -> reference.Hemisphere.proj(photo_stim_hemisphere="hemisphere")
    -> reference.CoordinateReference.proj(photo_stim_coordinate_ref="coordinate_ref")
    photo_stim_coordinate_ap: float    # in mm, anterior positive, posterior negative 
    photo_stim_coordinate_ml: float    # in mm, always postive, number larger when more lateral
    photo_stim_coordinate_dv: float    # in mm, always postive, number larger when more ventral (deeper)
    """

@schema
class Session(dj.Manual):
    definition = """
    -> subject.Subject
    session_time: datetime    # session time
    ---
    session_directory: varchar(256)
    session_note: varchar(256) # 
    """

    class Experimenter(dj.Part):
        definition = """
        -> master
        -> reference.Experimenter
        """

    class ExperimentType(dj.Part):
        definition = """
        -> master
        -> ExperimentType
        """

@schema
class Acquisition(dj.Manual):
    definition = """
    -> Session
    ---
    
    """
    class Behavior(dj.Part):
        definition = """
        -> master
        -> reference.BehavioralType
        ---
        behavior_sampling_rate: int
        behavior_time_stamp: longblob
        behavior_timeseries: longblob        
        """    
        
    class Ephys(dj.Part):
        definition = """
        -> master
        -> reference.EphysType
        ---
        ephys_sampling_rate: int
        ephys_time_stamp: longblob
        ephys_timeseries: longblob        
        """      
        
        
#@schema
#class TrialAcquisition(dj.Computed):
#    definition = """
#    -> Acquisition
#    -> behavior.TrialSet.Trial
#    ---
#    
#    """
#    class TrialBehavior(dj.Part):
#        definition = """
#
#        """
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    


























