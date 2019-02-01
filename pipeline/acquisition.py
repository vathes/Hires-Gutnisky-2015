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
        sess_data_dir = os.path.join('data', 'datafiles')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Trial import failed for session: {key["session_id"]}')
            return
        sess_data = sio.loadmat(os.path.join(sess_data_dir, 'data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat'),
                               struct_as_record = False, squeeze_me = True)['c']
        time_conversion_factor = utilities.time_unit_conversion_factor[
            sess_data.timeUnitNames[sess_data.timeSeriesArrayHash.value[0].timeUnit - 1]]  # (-1) to take into account Matlab's 1-based indexing
        behavior_data = sess_data.timeSeriesArrayHash.value[0].valueMatrix * time_conversion_factor
        time_stamps = sess_data.timeSeriesArrayHash.value[0].time * time_conversion_factor

        key['behavior_start_time'] = time_stamps[0]
        key['behavior_sampling_rate'] = 1/np.mean(np.diff(time_stamps))

        behavioral_key_maps = {'thetaAtBase': 'theta_at_base',
                               'amplitude': 'ampitude',
                               'phase': 'phase',
                               'setpoint': 'set_point',
                               'thetafilt': 'theta_filt',
                               'deltaKappa': 'delta_kappa',
                               'touch_onset': 'touch_onset',
                               'touch_offset': 'touch_offset',
                               'distance_to_pole': 'distance_to_pole',
                               'pole_available': 'pole_available',
                               'beam_break_times': 'beam_break_times'}
        self.insert1({**key, **{behavioral_key_maps[n]: behavior_data[n_idx, :]
                                for n_idx, n in enumerate(sess_data.timeSeriesArrayHash.value[0].idStr)}})
        print(f'Inserted behavioral data for session: {key["session_id"]}')


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
class EphysAcquisition(dj.Imported):
    definition = """ # Membrane potential recording from a cell, and electrical stimulation profile to this cell
    -> Cell
    ---
    voltage: longblob  # (mV)
    spike_train: longblob
    ephys_start_time: float  # (s)
    ephys_sampling_rate: float # (Hz)
    """     

    def make(self, key):
        sess_data_dir = os.path.join('data', 'datafiles')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Trial import failed for session: {key["session_id"]}')
            return
        sess_data = sio.loadmat(os.path.join(sess_data_dir, 'data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat'),
                               struct_as_record = False, squeeze_me = True)['c']
        time_conversion_factor = utilities.time_unit_conversion_factor[
            sess_data.timeUnitNames[sess_data.timeSeriesArrayHash.value[1].timeUnit - 1]]  # (-1) to take into account Matlab's 1-based indexing
        ephys_data = sess_data.timeSeriesArrayHash.value[1].valueMatrix * time_conversion_factor
        time_stamps = sess_data.timeSeriesArrayHash.value[1].time * time_conversion_factor

        key['voltage'] = ephys_data[0, :].todense()
        key['spike_train'] = ephys_data[1, :].todense()
        key['ephys_start_time'] = time_stamps[0]
        key['ephys_sampling_rate'] = 1/np.mean(np.diff(time_stamps))

        self.insert1(key)
        print(f'Inserted ephys data for session: {key["session_id"]}')


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
        sess_data_dir = os.path.join('data', 'datafiles')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Trial import failed for session: {key["session_id"]}')
            return
        sess_data = sio.loadmat(os.path.join(sess_data_dir, 'data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat'),
                               struct_as_record = False, squeeze_me = True)['c']
        key['trial_counts'] = len(sess_data.trialIds)
        self.insert1(key)

        # read trial info
        time_conversion_factor = utilities.time_unit_conversion_factor[
            sess_data.timeUnitNames[sess_data.trialTimeUnit - 1]]  # (-1) to take into account Matlab's 1-based indexing
        pole_in_times = sess_data.trialPropertiesHash.value[1] * time_conversion_factor
        pole_out_times = sess_data.trialPropertiesHash.value[2] * time_conversion_factor
        time_stamps = sess_data.timeSeriesArrayHash.value[1].time * time_conversion_factor

        print(f'Inserted trial set for session: {key["session_id"]}')
        print('Inserting trial ID: ', end="")
        for idx, trial_id in enumerate(sess_data.trialIds):
            key['trial_id'] = int(trial_id)
            key['start_time'] = time_stamps[np.where(
                sess_data.timeSeriesArrayHash.value[1].trial == int(trial_id))[0][0]]
            key['stop_time'] = time_stamps[np.where(
                sess_data.timeSeriesArrayHash.value[1].trial == int(trial_id))[0][-1]]
            key['trial_response'] = sess_data.trialTypeStr[np.where(sess_data.trialTypeMat[:-1, idx] == 1)[0][0]]
            key['trial_stim_present'] = int(sess_data.trialTypeMat[-1, idx] == 1)  # why DJ throws int type error for bool??
            key['trial_type'] = sess_data.trialPropertiesHash.value[-1][idx]
            key['pole_position'] = sess_data.trialPropertiesHash.value[0][idx] * 0.0992  # convert to micron here (0.0992 microns / microstep)
            self.Trial.insert1(key, ignore_extra_fields=True)
            # ======== Now add trial event timing to the EventTime part table ====
            event_dict = dict(trial_start=key['start_time'],
                              trial_stop=key['start_time'],
                              pole_in=pole_in_times[idx],
                              pole_out=pole_out_times[idx])
            self.EventTime.insert((dict(key, trial_event=k, event_time=v)
                                       for k, v in event_dict.items()),
                                      ignore_extra_fields=True)
            print(f'{trial_id} ', end = "")
        print('')
