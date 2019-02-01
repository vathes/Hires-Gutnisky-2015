'''
Schema of aquisition information.
'''
import re
import os
from datetime import datetime

import numpy as np
import scipy.io as sio
import datajoint as dj
import h5py as h5

from . import reference, subject, utilities, stimulation

schema = dj.schema(dj.config.get('database.prefix', '') + 'hg2015_acquisition')


@schema
class ExperimentType(dj.Lookup):
    definition = """
    experiment_type: varchar(64)
    """
    contents = [['behavior'], ['extracelluar'], ['photostim']]


@schema
class Session(dj.Manual):
    definition = """
    session_id: varchar(32)
    ---
    -> subject.Subject
    session_time: datetime    # session time
    session_directory = "": varchar(256)
    session_note = "" : varchar(256) 
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
class BehaviorAcquisition(dj.Imported):
    definition = """
    -> Session
    ---
    theta_at_base=null: longblob  #
    ampitude=null: longblob  #
    phase=null: longblob  #
    set_point=null: longblob  #
    theta_filt=null: longblob  #
    delta_kappa=null: longblob  #
    touch_onset=null: longblob  #
    touch_offset=null: longblob  #
    distance_to_pole=null: longblob  #
    pole_available=null: longblob  #
    beam_break_times=null: longblob  #
    behavior_start_time=null: float  #
    behavior_sampling_rate=null: float  #
    """

    def make(self, key):
        return


@schema
class PhotoStimulation(dj.Manual):
    definition = """ # Photostimulus profile used for stimulation in this session
    -> Session
    photostim_datetime: datetime # the time of performing this stimulation with respect to start time of the session, in the scenario of multiple stimulations per session
    ---
    -> stimulation.PhotoStimulationInfo
    photostim_timeseries=null: longblob  # (mW)
    photostim_start_time=null: float  # (s) first timepoint of photostim recording
    photostim_sampling_rate=null: float  # (Hz) sampling rate of photostim recording
    """    


@schema
class Cell(dj.Manual):
    definition = """ # A cell undergone intracellular recording in this session
    -> Session
    cell_id: varchar(36) # a string identifying the cell in which this intracellular recording is concerning
    ---
    cell_type: enum('excitatory', 'inhibitory', 'N/A')
    -> reference.BrainLocation
    recording_depth: decimal(6,2)  # (um)
    -> reference.WholeCellDevice
    """    
  
    
@schema
class IntracellularAcquisition(dj.Imported):
    definition = """ # Membrane potential recording from a cell, and electrical stimulation profile to this cell
    -> Cell
    """     
    
    class MembranePotential(dj.Part):
        definition = """
        -> master
        ---
        membrane_potential: longblob    # (mV) membrane potential recording at this cell
        membrane_potential_wo_spike: longblob # (mV) membrane potential without spike data, derived from membrane potential recording    
        membrane_potential_start_time: float # (s) first timepoint of membrane potential recording
        membrane_potential_sampling_rate: float # (Hz) sampling rate of membrane potential recording
        """
        
    class CurrentInjection(dj.Part):
        definition = """
        -> master
        ---
        current_injection: longblob
        current_injection_start_time: float  # first timepoint of current injection recording
        current_injection_sampling_rate: float  # (Hz) sampling rate of current injection recording
        """
        
    def make(self, key):
        return

    
@schema
class ProbeInsertion(dj.Manual):
    definition = """ # Description of probe insertion details during extracellular recording
    -> Session
    -> reference.Probe
    -> reference.ActionLocation
    """    
    

@schema
class ExtracellularAcquisition(dj.Imported):
    definition = """ # Raw extracellular recording, channel x time (e.g. LFP)
    -> ProbeInsertion
    """    
    
    class Voltage(dj.Part):
        definition = """
        -> master
        ---
        voltage: longblob   # (mV)
        voltage_start_time: float # (second) first timepoint of voltage recording
        voltage_sampling_rate: float # (Hz) sampling rate of voltage recording
        """
        
    def make(self, key):
        # this function implements the ingestion of raw extracellular data into the pipeline
        return None


@schema
class UnitSpikeTimes(dj.Imported):
    definition = """ 
    -> ProbeInsertion
    unit_id : smallint
    ---
    -> reference.Probe.Channel
    spike_times: longblob  # (s) time of each spike, with respect to the start of session 
    unit_cell_type: varchar(32)  # e.g. cell-type of this unit (e.g. wide width, narrow width spiking)
    unit_x: float  # (mm)
    unit_y: float  # (mm)
    unit_z: float  # (mm)
    spike_waveform: longblob  # waveform(s) of each spike at each spike time (spike_time x waveform_timestamps)
    """
        
    def make(self, key):
        # ================ Dataset ================
        sess_data_dir = os.path.join('..', 'data', 'extracellular', 'datafiles')
        # Get the Session definition from the keys of this session
        animal_id = key['subject_id']
        date_of_experiment = key['session_time']
        # Search the files in filenames to find a match for "this" session (based on key)
        sess_data_file = utilities.find_session_matched_nwbfile(sess_data_dir, animal_id, date_of_experiment)
        if sess_data_file is None: 
            print(f'UnitSpikeTimes import failed for: {animal_id} - {date_of_experiment}')
            return
        nwb = h5.File(os.path.join(sess_data_dir,sess_data_file), 'r')
        # ------ Spike ------
        ec_event_waveform = nwb['processing']['extracellular_units']['EventWaveform']
        ec_unit_times = nwb['processing']['extracellular_units']['UnitTimes']
        # - unit cell type
        cell_type = {}
        for tmp_str in ec_unit_times.get('cell_types').value:
            tmp_str = tmp_str.decode('UTF-8')
            split_str = re.split(' - ', tmp_str)
            cell_type[split_str[0]] = split_str[1]
        # - unit info
        print('Inserting spike unit: ', end="")
        for unit_str in ec_event_waveform.keys():
            unit_id = int(re.search('\d+', unit_str).group())
            unit_depth = ec_unit_times.get(unit_str).get('depth').value
            key['unit_id'] = unit_id
            key['channel_id'] = ec_event_waveform.get(unit_str).get('electrode_idx').value.item(0) - 1  # TODO: check if electrode_idx has MATLAB 1-based indexing (starts at 1)
            key['spike_times'] = ec_unit_times.get(unit_str).get('times').value
            key['unit_cell_type'] = cell_type[unit_str]
            key.update(zip(('unit_x', 'unit_y', 'unit_z'), unit_depth))
            key['spike_waveform'] = ec_event_waveform.get(unit_str).get('data').value
            self.insert1(key)
            print(f'{unit_id} ', end="")
        print('')
        nwb.close()
    

@schema
class TrialSet(dj.Imported):
    definition = """
    -> Session
    ---
    trial_counts: int # total number of trials
    """
    
    class Trial(dj.Part):
        definition = """
        -> master
        trial_id: smallint           # id of this trial in this trial set
        ---
        start_time = null: float               # start time of this trial, with respect to starting point of this session
        stop_time = null: float                # end time of this trial, with respect to starting point of this session
        -> reference.TrialType
        -> reference.TrialResponse
        trial_stim_present: bool  # is this a stim or no-stim trial
        pole_position: float  # (um)  the location of the pole along the anteroposterior axis of the animal
        """
        
    class EventTime(dj.Part):
        definition = """ # experimental paradigm event timing marker(s) for this trial
        -> master.Trial
        -> reference.ExperimentalEvent.proj(trial_event="event")
        ---
        event_time = null: float   # (in second) event time with respect to this session's start time
        """

    def make(self, key):
        sess_data_dir = os.path.join('..', 'data', 'datafiles')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Trial import failed for session: {key["session_id"]}')
            return
        sess_data = sio.loadmat(os.path.join(sess_data_dir, 'data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat'),
                               struct_as_record = False, squeeze_me = True)['c']
        # get time unit dict
        trial_time_unit = sess_data.timeUnitNames[sess_data.trialTimeUnit]




@schema
class TrialStimInfo(dj.Imported):
    definition = """ # information related to the stimulation settings for this trial
    -> TrialSet.Trial
    ---
    photo_stim_type: enum('stimulation','inhibition','N/A')
    photo_stim_period: enum('sample','delay','response','N/A')
    photo_stim_power: float  # (mW) stimulation power 
    photo_loc_galvo_x: float  # (mm) photostim coordinates field 
    photo_loc_galvo_y: float  # (mm) photostim coordinates field 
    photo_loc_galvo_z: float  # (mm) photostim coordinates field 
    """    
    
    def make(self, key):
        # this function implements the ingestion of Trial stim info into the pipeline
        return None
    
    