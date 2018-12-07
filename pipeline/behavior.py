'''
Schema of session information.
'''
import datajoint as dj
from pipeline import acquisition

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

        sessdatafile = 'data_structure_Cell'    
        
        print(sessdatafile)
        matfile = sio.loadmat(sessdatadir+sessdatafile, struct_as_record=False)
        sessdata = matfile['c'][0,0]
        header = re.split('_|\.',sessdatafile)
        
        animal_ID = header[3]
        date_of_experiment = header[4]
        cell = header[5]+'_'+header[6]
        
        animalId = Get1FromNestedArray(sessdata.animalId)
        sessdate = Get1FromNestedArray(sessdata.date)
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
        
        
        # ------------ TrialSet ------------
        if date_of_experiment is not None : 
            try: date_of_experiment = datetime.strptime(str(date_of_experiment),_datetimeformat_ymd) # expected datetime format: yymmdd
            except:
                try: date_of_experiment = datetime.strptime(str(date_of_experiment),_datetimeformat_ydm) # in case some dataset has messed up format: yyddmm
                except: date_of_experiment = None, print('Session Date error at: ' + date_of_experiment) # let's hope this doesn't happen
        
        
        TrialSet = behavior.TrialSet()
        
        behavior.TrialSet.connection.start_transaction()
        
        behavior.TrialSet.insert1(            
                            {'subject_id':animal_ID,
                             'cell_id':cell,
                             'session_time': date_of_experiment,
                             'trial_time_unit': timesUnitDict[trialTimeUnit],
                             'number_of_trials':len(trialIds)
                             }, 
                             skip_duplicates=True)
        
        for idx, trialId in enumerate(trialIds):
            tType = trialTypeMat[:,idx]
            tType = np.where(tType == 1)[0]
            tType = trialTypeStr[tType.item(0)] # this relies on the metadata consistency, e.g. a trial belongs to only 1 category of trial type
        
            this_trial_sample_idx = np.where(behav.trial[0,:] == (idx+1))[0] #  (+1) to take in to account that the native data format is MATLAB (index starts at 1)
            
            
            behavior.TrialSet.Trial.insert1(            
                            {'subject_id':animal_ID,
                             'cell_id':cell,
                             'session_time': date_of_experiment,
                             'trial_idx':trialId,
                             'trial_type':tType,
                             'pole_trial_condition': poleTrialCondition[idx],
                             'pole_position': polePos[idx],
                             'pole_in_time': poleInTime[idx],
                             'pole_out_time': poleOutTime[idx],
                             'lick_time': lickTime[0,idx],
                             'start_sample': this_trial_sample_idx[0],
                             'end_sample': this_trial_sample_idx[-1]
                             }, 
                             skip_duplicates=True)
        
        
        
        
        


















    


   
    
    