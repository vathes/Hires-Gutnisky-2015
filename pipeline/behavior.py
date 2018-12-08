'''
Schema of session information.
'''
import re
import os
from datetime import datetime
import numpy as np
import scipy.io as sio
import datajoint as dj
from pipeline import acquisition
from pipeline.helper_functions import Get1FromNestedArray, GetListFromNestedArray, _datetimeformat_ydm, _datetimeformat_ymd

schema = dj.schema('ttngu207_behavior',locals())

        
@schema
class TrialType(dj.Lookup):
    definition = """
    trial_type: varchar(64)
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
class TrialSet(dj.Imported):
    definition = """
    -> acquisition.Session
    ---
    number_of_trials: int
    trial_time_unit: enum('millisecond','second','minute','hour','day')
    """
    class Trial(dj.Part):
        definition = """
        -> master
        trial_idx: int
        ---
        -> TrialType
        pole_trial_condition: enum('Go','NoGo')  # string indicating whether the pole was presented in a ‘Go’ or ‘Nogo’ location
        pole_position: float                     # the location of the pole along the anteroposterior axis of the animal in microsteps (0.0992 microns / microstep)
        pole_in_time: float                      # the start of sample period for each trial (e.g. the onset of pole motion towards the exploration area), in units of seconds, relative to trialStartTimes
        pole_out_time: float                     # the end of the sample period (e.g. the onset of pole motion away from the exploration area)
        lick_time: longblob                      # an array of times of when the mouse’s tongue initiates contact with the spout
        
        start_sample: int       # the index of the starting sample of this trial, with respect to the starting of this session (at 0th sample point) - this way, time will be derived from the sampling rate of a particular recording downstream
        end_sample: int         # the index of the ending sample of this trial, with respect to the starting of this session (at 0th sample point) - this way, time will be derived from the sampling rate of a particular recording downstream
        
        """

    def _make_tuples(self,key):
        
        datadir = 'C://Users//thinh//Documents//TN-Vathes//NWB_Janelia_datasets//crcns_ssc5_data_HiresGutnisky2015//'
        sessdatadir = datadir + 'datafiles//'
        sessdatafiles = os.listdir(sessdatadir)
                
        # Get the Session definition from keys
        animal_id = key['subject_id']
        cell = key['cell_id']
        date_of_experiment = key['session_time']
                
        # Convert datetime to string format 
        date_of_experiment = datetime.strftime(date_of_experiment,_datetimeformat_ymd) # expected datetime format: yymmdd
        
        # Search the filenames to find a match for "this" session (based on key)
        sessdatafile = None
        for s in sessdatafiles:
            m1 = re.search(animal_id, s) 
            m2 = re.search(cell, s) 
            m3 = re.search(date_of_experiment, s) 
            if (m1 is not None) & (m2 is not None) & (m3 is not None):
                sessdatafile = s
                break
        
        # If session not found from dataset, break
        if sessdatafile is None:
            print(f'Session not found! - Subject: {animal_id} - Cell: {cell} - Date: {date_of_experiment}')
            return
        else: print(f'Found datafile: {sessdatafile}')
        
        # Now read the data and start ingesting
        matfile = sio.loadmat(sessdatadir+sessdatafile, struct_as_record=False)
        sessdata = matfile['c'][0,0]
        
        timeUnitIds = GetListFromNestedArray(sessdata.timeUnitIds)
        timeUnitNames = GetListFromNestedArray(sessdata.timeUnitNames)
        timesUnitDict = {}
        for idx, val in enumerate(timeUnitIds):
            timesUnitDict[val] = timeUnitNames[idx]
        
        trialIds = GetListFromNestedArray(sessdata.trialIds)
        trialStartTimes = GetListFromNestedArray(sessdata.trialStartTimes)
        trialTimeUnit = Get1FromNestedArray(sessdata.trialTimeUnit)
        trialTypeMat = sessdata.trialTypeMat
        trialTypeStr = GetListFromNestedArray(sessdata.trialTypeStr)
        
        trialPropertiesHash = sessdata.trialPropertiesHash[0,0]
        descr = GetListFromNestedArray(trialPropertiesHash.descr)
        keyNames = GetListFromNestedArray(trialPropertiesHash.keyNames)
        value = trialPropertiesHash.value
        polePos = GetListFromNestedArray(value[0,0])
        poleInTime = GetListFromNestedArray(value[0,1])
        poleOutTime = GetListFromNestedArray(value[0,2])
        lickTime = value[0,3]
        poleTrialCondition = GetListFromNestedArray(value[0,4])
        
        timeSeries = sessdata.timeSeriesArrayHash[0,0]
        behav = timeSeries.value[0,0][0,0]
        ephys = timeSeries.value[0,1][0,0]
        
                
        part_key = key.copy() # this is to perserve the original key for use in the part table later
        # form new key-values pair and insert key
        key['trial_time_unit'] = timesUnitDict[trialTimeUnit]
        key['number_of_trials'] = len(trialIds)
        self.insert1(key)
        print(f'Inserted trial set for session: Subject: {animal_id} - Cell: {cell} - Date: {date_of_experiment}')
        
        for idx, trialId in enumerate(trialIds):
            
            ### Debug here
            tmp = behav.trial[0,:]
            print('---')
            print(trialId)
            print(tmp)            
            ###
            
            tType = trialTypeMat[:,idx]
            tType = np.where(tType == 1)[0]
            tType = trialTypeStr[tType.item(0)] # this relies on the metadata consistency, e.g. a trial belongs to only 1 category of trial type
            
            try:
                this_trial_sample_idx = np.where(behav.trial[0,:] == trialId)[0] #  
            except: # this implementation is a safeguard against inconsistency in data formatting - e.g. "data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat" where "trial" vector is not referencing trialId
                this_trial_sample_idx = np.where(behav.trial[0,:] == (idx+1))[0] #  (+1) to take in to account that the native data format is MATLAB (index starts at 1)
            
            # form new key-values pair for part_key and insert
            part_key['trial_idx'] = trialId
            part_key['trial_type'] = tType
            part_key['pole_trial_condition'] = poleTrialCondition[idx]
            part_key['pole_position'] = polePos[idx]
            part_key['pole_in_time'] = poleInTime[idx]
            part_key['pole_out_time'] = poleOutTime[idx]
            part_key['lick_time'] = lickTime[0,idx]
            part_key['start_sample'] = this_trial_sample_idx[0]
            part_key['end_sample'] = this_trial_sample_idx[-1]
            self.Trial.insert1(part_key)
            print(f'Inserted trial ID: {trialId}')
            
        



        
        


















    


   
    
    