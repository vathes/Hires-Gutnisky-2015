import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

############## Explore #################
matdata = sio.loadmat('data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat', struct_as_record=False)

c = matdata['c'][0,0]
timeSeries = c.timeSeriesArrayHash[0,0]
behav = timeSeries.value[0,0][0,0]
ephys = timeSeries.value[0,1][0,0]

trialProperties = c.trialPropertiesHash[0,0]


#################################
import datajoint as dj
dj.config['database.host'] = 'tutorial-db.datajoint.io'

from pipeline import reference
from pipeline import subject
from pipeline import acquisition
from pipeline import behavior
from pipeline import ephys
from pipeline import action

all_erd = dj.ERD(reference) + dj.ERD(subject) + dj.ERD(action) + dj.ERD(acquisition) + dj.ERD(behavior) + dj.ERD(ephys)
all_erd.save('./all_erd.png')

############## INGESTION #################

# ==== Part 1: metadata ========

matdata = sio.loadmat('data_structure_Cell01_ANM244028_141021_JY1243_AAAA.mat', struct_as_record=False)

c = matdata['c'][0,0]
timeSeries = c.timeSeriesArrayHash[0,0]
behav = timeSeries.value[0,0][0,0]
ephys = timeSeries.value[0,1][0,0]

trialProperties = c.trialPropertiesHash[0,0]
    




















