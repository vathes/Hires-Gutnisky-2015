from pipeline import (reference, subject, acquisition, stimulation, analysis, virus,
                      intracellular, extracellular, behavior, utilities)

# ====================== Settings ======================
prioritized_sessions = acquisition.Session()
settings = {'reserve_jobs': True, 'suppress_errors': True, 'display_progress': False}

# ====================== Starting import and compute procedure ======================
# -- TrialSet
acquisition.TrialSet.populate(prioritized_sessions, **settings)
# -- Ephys
intracellular.MembranePotential.populate(prioritized_sessions, **settings)
intracellular.SpikeTrain.populate(prioritized_sessions, **settings)
# -- Behavioral
behavior.Behavior.populate(prioritized_sessions, **settings)
# -- Perform trial segmentation
print('------- Perform trial segmentation -------')
analysis.RealignedEvent.populate(**settings)
intracellular.TrialSegmentedMembranePotential.populate(prioritized_sessions, **settings)
intracellular.TrialSegmentedSpikeTrain.populate(prioritized_sessions, **settings)
behavior.TrialSegmentedBehavior.populate(prioritized_sessions, **settings)
stimulation.TrialSegmentedPhotoStimulus.populate(prioritized_sessions, **settings)

