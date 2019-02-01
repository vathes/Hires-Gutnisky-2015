import os
from datetime import datetime

import numpy as np
import scipy.io as sio
import re
from decimal import Decimal

import datajoint as dj
from pipeline import reference, subject, acquisition, stimulation, virus #, behavior, ephys, action
from pipeline import utilities

# Merge all schema and generate the overall ERD (then save in "/images")
all_erd = dj.ERD(reference) + dj.ERD(subject) + dj.ERD(acquisition)
all_erd.save('./images/all_erd.png')

# ================== Dataset ==================
data_dir = os.path.join('data')
meta_data_dir = os.path.join(data_dir, 'metadata')
sess_data_dir = os.path.join(data_dir, 'datafiles')

# ===================== Part 1: meta_data ========
meta_data_files = os.listdir(meta_data_dir)
for meta_data_file in meta_data_files:
    print(f'-- Read {meta_data_file} --')
    meta_data = sio.loadmat(os.path.join(
        meta_data_dir, meta_data_file), struct_as_record = False, squeeze_me=True)['meta_data']

    # ==================== subject ====================
    subject_info = {c: meta_data.__getattribute__(c)
                    for c in ('animal_ID', 'sex', 'species')}
    if meta_data.animal_background.size != 0:
        subject_info['subject_description'] = meta_data.animal_background
    # force subject_id to be lower-case for consistency
    subject_info['subject_id'] = subject_info.pop('animal_ID').lower()
    # dob and sex
    subject_info['sex'] = subject_info['sex'][0].upper() if subject_info['sex'].size != 0 else 'U'
    if meta_data.data_of_birth.size != 0:
        subject_info['date_of_birth'] = utilities.parse_prefix(meta_data.data_of_birth)

    # source and strain
    source_strain = meta_data.source_strain  # subject.Allele
    if len(source_strain) > 0:  # if found, search found string to find matched strain in db
        for s in subject.StrainAlias.fetch():
            m = re.search(re.escape(s[0]), source_strain, re.I)
            if m is not None:
                subject_info['strain'] = (subject.StrainAlias & {'strain_alias': s[0]}).fetch1('strain')
                break
    else:
        subject_info['strain'] = 'N/A'

    source_identifier = meta_data.source_identifier  # reference.AnimalSource
    if len(source_identifier) > 0:  # if found, search found string to find matched strain in db
        for s in reference.AnimalSourceAlias.fetch():
            m = re.search(re.escape(s[0]), source_identifier, re.I)
            if m is not None:
                subject_info['animal_source'] = (reference.AnimalSourceAlias & {'animal_source_alias': s[0]}).fetch1(
                    'animal_source')
                break
    else:
        subject_info['animal_source'] = 'N/A'

    if subject_info not in subject.Subject.proj():
        subject.Subject.insert1(subject_info, ignore_extra_fields = True)

    # ==================== session ====================
    # -- session_time
    date_of_experiment = utilities.parse_prefix(str(meta_data.date_of_experiment))  # acquisition.Session
    if date_of_experiment is not None:
        session_info = {'session_id': re.search('Cell\d+', meta_data_file).group(),
                       'session_time': date_of_experiment}
        # experimenter and experiment type (possible multiple experimenters or types)
        experiment_types = meta_data.experiment_type
        experimenters = meta_data.experimenters  # reference.Experimenter
        experimenters = [experimenters] if np.array(experimenters).size <= 1 else experimenters  # in case there's only 1 experimenter

        reference.Experimenter.insert(({'experimenter': k} for k in experimenters
                                        if {'experimenter': k} not in reference.Experimenter))
        acquisition.ExperimentType.insert(({'experiment_type': k} for k in experiment_types
                                        if {'experiment_type': k} not in acquisition.ExperimentType))

        if {**subject_info, **session_info} not in acquisition.Session.proj():
            with acquisition.Session.connection.transaction:
                acquisition.Session.insert1({**subject_info, **session_info}, ignore_extra_fields=True)
                acquisition.Session.Experimenter.insert((dict({**subject_info, **session_info}, experimenter=k) for k in experimenters), ignore_extra_fields=True)
                acquisition.Session.ExperimentType.insert((dict({**subject_info, **session_info}, experiment_type=k) for k in experiment_types), ignore_extra_fields=True)
            print(f'Creating Session - Subject: {subject_info["subject_id"]} - Date: {session_info["session_time"]}')

    # ==================== Intracellular ====================
    if isinstance(meta_data.extracellular, sio.matlab.mio5_params.mat_struct):
        extracellular = meta_data.extracellular
        brain_region = re.split(',\s?|\s', extracellular.atlas_location)[1]
        recording_coord_depth = extracellular.recording_coord_location[1]  # acquisition.RecordingLocation
        cortical_layer, brain_subregion = re.split(',\s?|\s', extracellular.recording_coord_location[0])[1:3]
        hemisphere = 'left'  # hardcoded here, not found in data, not found in paper
        brain_location = {'brain_region': brain_region,
                          'brain_subregion': brain_subregion,
                          'cortical_layer': cortical_layer,
                          'hemisphere': hemisphere,
                          'brain_location_full_name': extracellular.atlas_location}
        # -- BrainLocation
        if brain_location not in reference.BrainLocation.proj():
            reference.BrainLocation.insert1(brain_location)
        # -- Whole Cell Device
        ie_device = extracellular.probe_type.split(', ')[0]
        if {'device_name': ie_device} not in reference.WholeCellDevice.proj():
            reference.WholeCellDevice.insert1({'device_name': ie_device, 'device_desc': extracellular.probe_type})
        # -- Cell
        cell_id = meta_data.cell.upper()
        cell_key = dict({**subject_info, **session_info, **brain_location},
                        cell_id=cell_id,
                        cell_type=extracellular.cell_type,
                        device_name=ie_device,
                        recording_depth=round(Decimal(re.match(
                            '\d+', extracellular.recording_coord_location[1]).group()), 2))
        if cell_key not in acquisition.Cell.proj():
            acquisition.Cell.insert1(cell_key, ignore_extra_fields = True)
            print(f'\tInsert Cell: {cell_id}')

    # ==================== Photo stimulation ====================
    if isinstance(meta_data.photostim, sio.matlab.mio5_params.mat_struct):
        photostimInfo = meta_data.photostim
        brain_region = re.split(',\s?|\s', photostimInfo.photostim_atlas_location)[1]
        coord_ap_ml_dv = re.findall('\d+.\d+', photostimInfo.photostim_coord_location)

        brain_location = {'brain_region': brain_region,
                          'brain_subregion': 'N/A',
                          'cortical_layer': 'N/A',
                          'hemisphere': hemisphere,
                          'brain_location_full_name': photostimInfo.photostim_atlas_location}
        # -- BrainLocation
        if brain_location not in reference.BrainLocation.proj():
            reference.BrainLocation.insert1(brain_location)
        # -- ActionLocation
        action_location = dict(brain_location,
                               coordinate_ref = 'bregma',
                               coordinate_ap = round(Decimal(coord_ap_ml_dv[0]), 2),
                               coordinate_ml = round(Decimal(coord_ap_ml_dv[1]), 2),
                               coordinate_dv = round(Decimal('0'), 2))  # no depth information for photostim
        if action_location not in reference.ActionLocation.proj():
            reference.ActionLocation.insert1(action_location, ignore_extra_fields=True)

        # -- Device
        stim_device = 'laser'  # hard-coded here..., could not find a more specific name from metadata
        if {'device_name': stim_device} not in stimulation.PhotoStimDevice.proj():
            stimulation.PhotoStimDevice.insert1({'device_name': stim_device})

        # -- PhotoStimulationInfo
        photim_stim_info = dict(action_location,
                                device_name=stim_device,
                                photo_stim_excitation_lambda=float(re.match(
                                    '\d+', photostimInfo.__getattribute__('lambda')).group()),
                                photo_stim_method=photostimInfo.stimulation_method)
        if photim_stim_info not in stimulation.PhotoStimulationInfo.proj():
            stimulation.PhotoStimulationInfo.insert1(photim_stim_info, ignore_extra_fields=True)
            print(f'\tCreate new Photostim Information')

        if dict(session_info, photostim_datetime = session_info['session_time']) not in acquisition.PhotoStimulation.proj():
            acquisition.PhotoStimulation.insert1(dict({**subject_info, **session_info, **photim_stim_info},
                                                      photostim_datetime = session_info['session_time']),
                                                 ignore_extra_fields = True)
            print(f'\tInsert Photostim')

    # ==================== Virus ====================
    if isinstance(meta_data.virus, sio.matlab.mio5_params.mat_struct):
        virus_info = dict(
            virus_source=meta_data.virus.virus_source,
            virus=meta_data.virus.virus_name,
            virus_lot_number=meta_data.virus.virus_lot_number if meta_data.virus.virus_lot_number.size != 0 else '',
            virus_titer=meta_data.virus.titer.replace('x10', '') if len(meta_data.virus.titer) > 0 else None)

        if virus_info not in virus.Virus.proj():
            virus.Virus.insert1(virus_info)

        brain_location = {'brain_region': meta_data.virus.atlas_location.split(' ')[0],
                          'brain_subregion': meta_data.virus.virus_coord_location,
                          'cortical_layer': 'N/A',
                          'hemisphere': hemisphere}
        # -- BrainLocation
        if brain_location not in reference.BrainLocation.proj():
            reference.BrainLocation.insert1(brain_location)

        virus_injection = dict(
            {**virus_info, **subject_info, **brain_location},
            coordinate_ref='bregma',
            injection_date=utilities.parse_prefix(meta_data.virus.injection_date))

        virus.VirusInjection.insert([dict(virus_injection,
                                          injection_depth = round(Decimal(re.match('\d+', depth).group()), 2),
                                          injection_volume = round(Decimal(re.match('\d+', vol).group()), 2))
                                     for depth, vol in zip(meta_data.virus.depth, meta_data.virus.volume)],
                                    ignore_extra_fields=True, skip_duplicates=True)
        print(f'\tInsert Virus Injections - Count: {len(meta_data.virus.depth)}')

# ====================== Starting import and compute procedure ======================
# -- TrialSet
acquisition.TrialSet.populate()
# -- Ephys
acquisition.EphysAcquisition.populate()
# -- Behavioral
acquisition.BehaviorAcquisition.populate()



# ========================================================================
# =========================== R&D ==========================
# sess_data = sio.loadmat(os.path.join(os.path.join('data', 'datafiles'), 'data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat'),
#                       struct_as_record=False, squeeze_me=True)['c']




# =========================== NOT INGESTED ==========================
# source_gene_copy = meta_data.source_gene_copy  # probably  reference.SourceStrain (not implemented)
# source_transgene = meta_data.source_transgene  # not sure
#
# fiber = meta_data.fiber  # empty
# manipulation_type = meta_data.manipulation_type  # not sure
# onsetLatency = meta_data.onsetLatency  # not sure
# weight_after_experiment = meta_data.weight_after_experiment  # action.Weighing
# weight_before_experiment = meta_data.weight_before_experiment  # action.Weighing
# whisker = meta_data.whisker  # reference.Whisker
#
# extracellular = meta_data.extracellular if 'extracellular' in meta_data._fieldnames else None
# ground_location = meta_data.extracellular.ground_location  # not sure
# identification_method = meta_data.extracellular.identification_method  # not sure
# nth_time_accessing_tissue = meta_data.extracellular.nth_time_accessing_tissue
# recording_location_marker = meta_data.extracellular.recording_location_marker  # not sure
# recording_type = meta_data.extracellular.recording_type  # not sure
# ref_loc = meta_data.extracellular.ref_loc  # not sure
#
# injection_pattern = meta_data.virus.injection_pattern  # empty         # action.VirusInjection
# task_keyword = meta_data.behaviorInfo.task_keyword