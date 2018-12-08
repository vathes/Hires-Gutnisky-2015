import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import os
from datetime import datetime
import re

from pipeline.helper_functions import Get1FromNestedArray, GetListFromNestedArray, _datetimeformat_ydm, _datetimeformat_ymd

############## Dataset #################

datadir = 'C://Users//thinh//Documents//TN-Vathes//NWB_Janelia_datasets//crcns_ssc5_data_HiresGutnisky2015//'
metadatadir = datadir + 'metadata//'
sessdatadir = datadir + 'datafiles//'

#################################
import datajoint as dj
dj.config['database.host'] = 'tutorial-db.datajoint.io'

import pipeline
from pipeline import reference, subject, acquisition, behavior, ephys, action

all_erd = dj.ERD(reference) + dj.ERD(subject) + dj.ERD(action) + dj.ERD(acquisition) + dj.ERD(behavior) + dj.ERD(ephys)
all_erd.save('./images/all_erd.png')

############## INGESTION #################

# ===================== Part 1: metadata ========

metadatafiles = os.listdir(metadatadir)
for metadatafile in metadatafiles:    
    print(metadatafile)
    matfile = sio.loadmat(metadatadir+metadatafile, struct_as_record=False)
    metadata = matfile['meta_data'][0,0]
    
    animal_ID = Get1FromNestedArray(metadata.animal_ID).upper()            # subject.Subject # .upper() here to handle some inconsistencies in the dataset
    animal_background = Get1FromNestedArray(metadata.animal_background)    # empty
    cell = Get1FromNestedArray(metadata.cell)                              # reference.Cell
    dob = Get1FromNestedArray(metadata.data_of_birth)                      # subject.Subject
    date_of_experiment = Get1FromNestedArray(metadata.date_of_experiment)  # acquisition.Session
    
    experiment_type = metadata.experiment_type
    
    experimenters = Get1FromNestedArray(metadata.experimenters)            # reference.Experimenter
    fiber = metadata.fiber # empty
    manipulation_type = Get1FromNestedArray(metadata.manipulation_type)    # not sure
    onsetLatency = Get1FromNestedArray(metadata.onsetLatency)              # not sure
    sex = Get1FromNestedArray(metadata.sex)                                # subject.Subject
    source_gene_copy = Get1FromNestedArray(metadata.source_gene_copy)      # probably  reference.SourceStrain (not implemented)
    source_identifier = Get1FromNestedArray(metadata.source_identifier)    # reference.AnimalSource
    source_strain = Get1FromNestedArray(metadata.source_strain)            # subject.Allele
    source_transgene = Get1FromNestedArray(metadata.source_transgene)      # not sure
    species = Get1FromNestedArray(metadata.species)                        # subject.Species
    weight_after_experiment = Get1FromNestedArray(metadata.weight_after_experiment)   # action.Weighing
    weight_before_experiment = Get1FromNestedArray(metadata.weight_before_experiment) # action.Weighing
    whisker = Get1FromNestedArray(metadata.whisker)                        # reference.Whisker
    
    extracellular = metadata.extracellular[0,0] if metadata.extracellular.size != 0 else None  
    if  extracellular is not None:
        atlas_location = Get1FromNestedArray(extracellular.atlas_location)       # not sure
        cell_type = Get1FromNestedArray(extracellular.cell_type)                 # reference.Cell
        ground_location =  Get1FromNestedArray(extracellular.ground_location)    # not sure 
        identification_method =  Get1FromNestedArray(extracellular.identification_method) # not sure
        nth_time_accessing_tissue =  Get1FromNestedArray(extracellular.nth_time_accessing_tissue)
        probe_type =  Get1FromNestedArray(extracellular.probe_type)              # probably reference.Device (not implemented)
        recording_coord_location = Get1FromNestedArray(extracellular.recording_coord_location[0,0]) # reference.BrainLocation
        recording_coord_depth =    Get1FromNestedArray(extracellular.recording_coord_location[0,1]) # acquisition.RecordingLocation 
        recording_location_marker =  Get1FromNestedArray(extracellular.recording_location_marker) # not sure
        recording_type = Get1FromNestedArray(extracellular.recording_type)       # not sure
        ref_loc =  Get1FromNestedArray(extracellular.ref_loc)                            # not sure
    else: 
        atlas_location, cell_type, ground_location, identification_method, nth_time_accessing_tissue, \
        probe_type, recording_coord_location, recording_coord_depth, recording_location_marker, \
        recording_type, ref_loc = None, None, None, None, None, None, None, None, None, None, None 
    
    behaviorInfo = metadata.behavior[0,0] if metadata.behavior.size != 0 else None  
    if  behaviorInfo is not None:
        task_keyword = behaviorInfo.task_keyword                        # not sure
    else: task_keyword = None
    
    photostim = metadata.photostim[0,0] if metadata.photostim.size != 0 else None    
    if  photostim is not None:
        #stim_lambda =                                                        # acquisition.PhotoStim
        photostim_atlas_location = Get1FromNestedArray(photostim.photostim_atlas_location) # acquisition.PhotoStim
        photostim_coord_location = Get1FromNestedArray(photostim.photostim_coord_location) # acquisition.PhotoStim
        stimulation_method = Get1FromNestedArray(photostim.stimulation_method)
    else: photostim_atlas_location, photostim_coord_location, stimulation_method = None, None, None
    
    virus = metadata.virus[0,0] if metadata.virus.size != 0 else None  
    if  virus is not None:
        atlas_location = Get1FromNestedArray(virus.atlas_location)               # not sure (location)
        depth = virus.depth                                         # action.VirusInjection
        injection_date = Get1FromNestedArray(virus.injection_date)               # action.VirusInjection
        injection_pattern = virus.injection_pattern # empty         # action.VirusInjection
        titer = Get1FromNestedArray(virus.titer)                                 # reference.Virus
        virus_coord_location = Get1FromNestedArray(virus.virus_coord_location)  # not sure (location)
        virus_lot_number = virus.virus_lot_number # empty           # reference.Virus
        virus_name = Get1FromNestedArray(virus.virus_name)                       # reference.Virus
        virus_source = Get1FromNestedArray(virus.virus_source)                   # reference.VirusSource
        volume = virus.volume                                       # action.VirusInjection
    else: 
        atlas_location, depth, injection_date, injection_pattern, titer, \
        virus_coord_location, virus_lot_number, virus_name, virus_source, \
        volume = None, None, None, None, None, None, None, None, None, None
        
    
    
    # ------------ Species ------------
    subject.Species.insert1([species], skip_duplicates=True)
    
    # ------------ Strain ------------
    subject.Strain.insert1(['N/A'], skip_duplicates=True)
    
    # ------------ Allele ------------
    subject.Allele().insert1([source_strain], skip_duplicates=True)
    
    # ------------ Animal Source ------------
    if source_identifier is None : source_identifier = 'N/A'
    reference.AnimalSource.insert1([source_identifier], skip_duplicates=True)
    
    # ------------ Subject ------------
    if sex is None : sex = 'U'
    if dob is not None : 
        try: dob = datetime.strptime(str(dob),_datetimeformat_ymd) # expected datetime format: yymmdd
        except:
            try: dob = datetime.strptime(str(dob),_datetimeformat_ydm) # in case some dataset has messed up format: yyddmm
            except: dob = None
            
    strain = 'N/A' # for this dataset, let's just say we dont know strain
    
    if dob is not None:
        subject.Subject.insert1(
                {'subject_id':animal_ID,
                 'species':species,
                 'strain': strain,
                 'animal_source':source_identifier,
                 'sex': sex,
                 'date_of_birth': dob}, 
                 skip_duplicates=True)
    else: 
        subject.Subject.insert1(
                {'subject_id':animal_ID,
                 'species':species,
                 'strain': strain,
                 'animal_source':source_identifier,
                 'sex': sex},
                 skip_duplicates=True)
    
    # ------------ Cell ------------
    
    # handling some inconsistency in cell naming convention in the metadata
    if len(cell) < 7 : cell = cell + '_AAAA'
    
    reference.Cell.insert1([cell, cell_type],skip_duplicates=True)
    
    # ------------ RecordingLocation ------------
    splittedstr = re.split(', |Layer | ',recording_coord_location)
    layer = splittedstr[1].lower() 
    subregion = splittedstr[2].lower()  
    region = splittedstr[3].lower()  
    depth = re.findall(r'\d+',recording_coord_depth)[0] # need to implement a way to check if the unit is in um, if not then convert
    
    # Check brain location first
    matchedBrainLocation = (reference.BrainLocation() & {'brain_location': region , 'cortical_layer' : layer , 'brain_subregion' : subregion}).fetch()
    if matchedBrainLocation.size == 0: 
        reference.BrainLocation.insert1(
                {'brain_location':region,
                 'cortical_layer':layer,
                 'brain_subregion':subregion},
                skip_duplicates=True)        
        
    
    acquisition.RecordingLocation.insert1(
            {'brain_location':region,
             'cortical_layer':layer,
             'brain_subregion':subregion,
             'recording_depth':depth},
            skip_duplicates=True)
    
    # ------------ Experimenter ------------
    reference.Experimenter.insert1([experimenters],skip_duplicates=True)
    
    # ------------ Session ------------
    if date_of_experiment is not None : 
        try: date_of_experiment = datetime.strptime(str(date_of_experiment),_datetimeformat_ymd) # expected datetime format: yymmdd
        except:
            try: date_of_experiment = datetime.strptime(str(date_of_experiment),_datetimeformat_ydm) # in case some dataset has messed up format: yyddmm
            except: date_of_experiment = None, print('Session Date error at: ' + metadatafile) # let's hope this doesn't happen
        
    
    if date_of_experiment is not None: 
        acquisition.Session.connection.start_transaction()
        acquisition.Session.insert1(            
                    {'subject_id':animal_ID,
                     'cell_id':cell,
                     'session_time': date_of_experiment,
                     'brain_location':region,
                     'cortical_layer':layer,
                     'brain_subregion':subregion,
                     'recording_depth':depth,
                     }, 
                     skip_duplicates=True)
        acquisition.Session.Experimenter.insert1(            
                    {'subject_id':animal_ID,
                     'cell_id':cell,
                     'session_time': date_of_experiment,
                     'experimenter': experimenters
                     }, 
                     skip_duplicates=True)
        # there is still the ExperimentType part table here...
        acquisition.Session.connection.commit_transaction()
        print(f'\tSession created - Subject: {animal_ID} - Cell: {cell} - Date: {date_of_experiment}')
    
    
    ## Need to perform ingestion for Virus, VirusInjection and PhotoStim ##
    


# ===================== Part 2: Acquisition data (still in experimental stage) ========

behavior.TrialSet.populate()



#=========================== R&D ==========================

matdata = sio.loadmat('data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat', struct_as_record=False)

c = matdata['c'][0,0]
timeSeries = c.timeSeriesArrayHash[0,0]
behav = timeSeries.value[0,0][0,0]
ephys = timeSeries.value[0,1][0,0]

trialProperties = c.trialPropertiesHash[0,0]




































