from pipeline import (reference, subject, acquisition, stimulation, analysis, virus,
                      intracellular, extracellular, behavior, utilities)

# ====================== Settings ======================
prioritized_sessions = acquisition.Session()

reserve_jobs = True
suppress_errors = True

# ====================== Starting import and compute procedure ======================
# -- TrialSet
acquisition.TrialSet.populate(prioritized_sessions,
                              reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
# -- Ephys
intracellular.MembranePotential.populate(prioritized_sessions,
                                         reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
intracellular.SpikeTrain.populate(prioritized_sessions,
                                  reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
# -- Behavioral
behavior.Behavior.populate(prioritized_sessions,
                           reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
# -- Perform trial segmentation
print('------- Perform trial segmentation -------')
analysis.RealignedEvent.populate(reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
intracellular.TrialSegmentedMembranePotential.populate(prioritized_sessions,
                                                       reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
intracellular.TrialSegmentedSpikeTrain.populate(prioritized_sessions,
                                                reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
behavior.TrialSegmentedBehavior.populate(prioritized_sessions,
                                         reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
stimulation.TrialSegmentedPhotoStimulus.populate(prioritized_sessions,
                                                 reserve_jobs=reserve_jobs, suppress_errors=suppress_errors)
