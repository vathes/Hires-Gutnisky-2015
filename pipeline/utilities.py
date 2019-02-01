import os
from datetime import datetime
import re

import h5py as h5
import numpy as np

from . import reference, acquisition


# datetime format - should probably read this from a config file and not hard coded here
datetimeformat_ymdhms = '%Y-%m-%d %H:%M:%S'
datetimeformat_ymd = '%y%m%d'
datetimeformat_ydm = '%y%d%m'

time_unit_conversion_factor = {
        'millisecond': 10e-3,
        'second': 1,
        'minute': 60,
        'hour': 3600,
        'day': 86400
        }

def parse_prefix(line):
    cover = len(datetime.now().strftime(datetimeformat_ymd))
    try:
        return datetime.strptime(line[:cover], datetimeformat_ymd)
    except Exception as e:
        msg = f'Error:  {str(e)} \n'
        cover = len(datetime.now().strftime(datetimeformat_ydm))
        try:
            return datetime.strptime(line[:cover], datetimeformat_ydm)
        except Exception as e:
            print(f'{msg}\t{str(e)}\n\tReturn None')
            return None    


def find_session_matched_matfile(sess_data_dir, sess_key):
        sess_data_files = os.listdir(sess_data_dir)
        # Search the filenames to find a match for "this" session (based on key)
        sess_data_file = None
        for s in sess_data_files:
            if re.search(sess_key['session_id'], s):
                sess_data_file = s
                print(f'Found datafile: {sess_data_file}')
        if sess_data_file:
            return sess_data_file
        else:
            print(f'Session not found! - Session: {sess_key["session_id"]}')
            return None
        