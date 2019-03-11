#!/usr/bin/env python3
import os

import sys
from datetime import datetime
from dateutil.tz import tzlocal
import pytz
import re
import numpy as np
import pandas as pd

from pipeline import (reference, subject, acquisition, stimulation, analysis, virus,
                      intracellular, extracellular, behavior, utilities)
import pynwb
from pynwb import NWBFile, NWBHDF5IO

# =============================================
# Each NWBFile represent a session, thus for every session in acquisition.Session, we build one NWBFile

for session_key in acquisition.Session.fetch('KEY'):
    this_session = (acquisition.Session & session_key).fetch1()
    # =============== General ====================
    # -- NWB file - a NWB2.0 file for each session
    nwbfile = NWBFile(
        session_description=this_session['session_note'],
        identifier='_'.join([this_session['subject_id'],
                             this_session['session_time'].strftime('%Y-%m-%d_%H-%M-%S')]),
        session_start_time=this_session['session_time'],
        file_create_date=datetime.now(tzlocal()),
        experimenter='; '.join((acquisition.Session.Experimenter & session_key).fetch('experimenter')),
        institution='Janelia Research Campus',  # TODO: not in pipeline
        related_publications='doi:10.1038/nature22324')  # TODO: not in pipeline
    # -- subject
    subj = (subject.Subject & session_key).fetch1()
    nwbfile.subject = pynwb.file.Subject(
        subject_id=this_session['subject_id'],
        description=subj['subject_description'],
        genotype= ' x '.join((subject.Subject & session_key).fetch('allele')),
        sex=subj['sex'],
        species=subj['species'])
    # =============== Intracellular ====================
    cell = ((acquisition.Cell & session_key).fetch1()
            if len(acquisition.Cell & session_key) == 1
            else None)
    if cell:
        # metadata
        whole_cell_device = nwbfile.create_device(name=cell['device_name'])
        ic_electrode = nwbfile.create_ic_electrode(
            name=cell['cell_id'],
            device=whole_cell_device,
            description='N/A',
            filtering='low-pass: 10kHz',  # TODO: not in pipeline
            location='; '.join([f'{k}: {str(v)}'
                                for k, v in (reference.ActionLocation & cell).fetch1().items()]))
        # acquisition - membrane potential
        mp, mp_timestamps = (intracellular.MembranePotential & cell).fetch1(
            'membrane_potential', 'membrane_potential_timestamps')
        nwbfile.add_acquisition(pynwb.icephys.PatchClampSeries(name='membrane_potential',
                                                               electrode=ic_electrode,
                                                               unit='mV',  # TODO: not in pipeline
                                                               conversion=1e-3,
                                                               gain=1.0,  # TODO: not in pipeline
                                                               data=mp,
                                                               timestamps=mp_timestamps))

        # acquisition - spike train
        spk, spk_timestamps = (intracellular.SpikeTrain & cell).fetch1(
            'spike_train', 'spike_timestamps')
        nwbfile.add_acquisition(pynwb.icephys.PatchClampSeries(name='spike_train',
                                                               electrode=ic_electrode,
                                                               unit='a.u.',  # TODO: not in pipeline
                                                               conversion=1,
                                                               gain=1.0,
                                                               data=spk,
                                                               timestamps=spk_timestamps))

    # =============== Behavior ====================
    behavior_data = ((behavior.Behavior & session_key).fetch1()
                       if len(behavior.Behavior & session_key) == 1
                       else None)
    if behavior_data:
        behav_acq = pynwb.behavior.BehavioralTimeSeries(name = 'behavior')
        nwbfile.add_acquisition(behav_acq)
        [behavior_data.pop(k) for k in behavior.Behavior.primary_key]
        timestamps = behavior_data.pop('behavior_timestamps')
        for b_k, b_v in behavior_data.items():
            behav_acq.create_timeseries(name=b_k,
                                        unit='a.u.',
                                        conversion=1.0,
                                        data=b_v,
                                        timestamps=timestamps)

    # =============== Photostimulation ====================
    photostim = ((acquisition.PhotoStimulation & session_key).fetch1()
                       if len(acquisition.PhotoStimulation & session_key) == 1
                       else None)
    if photostim:
        photostim_device = (stimulation.PhotoStimDevice & photostim).fetch1()
        stim_device = nwbfile.create_device(name=photostim_device['device_name'])
        stim_site = pynwb.ogen.OptogeneticStimulusSite(
            name='-'.join([photostim['hemisphere'], photostim['brain_region']]),
            device=stim_device,
            excitation_lambda=float(photostim['photo_stim_excitation_lambda']),
            location = '; '.join([f'{k}: {str(v)}' for k, v in
                                  (reference.ActionLocation & photostim).fetch1().items()]),
            description=(stimulation.PhotoStimulationInfo & photostim).fetch1('photo_stim_notes'))
        nwbfile.add_ogen_site(stim_site)

        if photostim['photostim_timeseries'] is not None:
            nwbfile.add_stimulus(pynwb.ogen.OptogeneticSeries(
                name='_'.join(['photostim_on', photostim['photostim_datetime'].strftime('%Y-%m-%d_%H-%M-%S')]),
                site=stim_site,
                unit = 'mW',
                resolution = 0.0,
                conversion = 1e-6,
                data = photostim['photostim_timeseries'],
                starting_time = photostim['photostim_start_time'],
                rate = photostim['photostim_sampling_rate']))

    # =============== TrialSet ====================
    # NWB 'trial' (of type dynamic table) by default comes with three mandatory attributes:
    #                                                                       'id', 'start_time' and 'stop_time'.
    # Other trial-related information needs to be added in to the trial-table as additional columns (with column name
    # and column description)
    if len(acquisition.TrialSet & session_key).fetch() == 1:
        # Get trial descriptors from TrialSet.Trial and TrialStimInfo
        trial_columns = [{'name': tag,
                          'description': re.sub('\s+:|\s+', ' ', re.search(
                              f'(?<={tag})(.*)', str((acquisition.TrialSet.Trial * stimulation.TrialPhotoStimInfo).heading)).group())}
                         for tag in (acquisition.TrialSet.Trial * stimulation.TrialPhotoStimInfo).fetch(as_dict=True, limit=1)[0].keys()
                         if tag not in (acquisition.TrialSet.Trial & stimulation.TrialPhotoStimInfo).primary_key + ['start_time', 'stop_time']]

        # Trial Events
        trial_events = set((acquisition.TrialSet.EventTime & session_key).fetch('trial_event'))
        event_names = [{'name': e, 'description': d}
                       for e, d in zip(*(reference.ExperimentalEvent & [{'event': k}
                                                                        for k in trial_events]).fetch('event',
                                                                                                      'description'))]
        # Add new table columns to nwb trial-table for trial-label
        for c in trial_columns + event_names:
            nwbfile.add_trial_column(**c)

        photostim_tag_default = {tag: '' for tag in stimulation.TrialPhotoStimInfo().fetch(as_dict=True, limit=1)[0].keys()
                                 if tag not in stimulation.TrialPhotoStimInfo.primary_key}
        # Add entry to the trial-table
        for trial in (acquisition.TrialSet.Trial & session_key).fetch(as_dict=True):
            events = dict(zip(*(acquisition.TrialSet.EventTime & trial).fetch('trial_event', 'event_time')))

            photostim_tag = (stimulation.TrialPhotoStimInfo & trial).fetch(as_dict=True)
            trial_tag_value = {**trial, **events, **photostim_tag[0]} if len(photostim_tag) == 1 else {**trial, **photostim_tag_default}
            # rename 'trial_id' to 'id'
            trial_tag_value['id'] = trial_tag_value['trial_id']
            [trial_tag_value.pop(k) for k in acquisition.TrialSet.Trial.primary_key]
            nwbfile.add_trial(**trial_tag_value)

    # =============== Write NWB 2.0 file ===============
    save_path = os.path.join('data', 'NWB 2.0')
    save_file_name = ''.join([nwbfile.identifier, '.nwb'])
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with NWBHDF5IO(os.path.join(save_path, save_file_name), mode = 'w') as io:
        io.write(nwbfile)
        print(f'Write NWB 2.0 file: {save_file_name}')






