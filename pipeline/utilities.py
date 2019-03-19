import os
from datetime import datetime
import re

import h5py as h5
import numpy as np

insert_size = 10

time_unit_conversion_factor = {'millisecond': 1e-3,
                               'second': 1,
                               'minute': 60,
                               'hour': 3600,
                               'day': 86400}

datetime_formats = ('%y%m%d', '%y%d%m')


def parse_date(text):
    for fmt in datetime_formats:
        cover = len(datetime.now().strftime(fmt))
        try:
            return datetime.strptime(text[:cover], fmt)
        except ValueError:
            pass
    raise ValueError('no valid date format found')


def find_session_matched_matfile(sess_data_dir, sess_key):
    """
    Search the filenames to find a match for "this" session (based on key)
    """
    sess_data_files = os.listdir(sess_data_dir)
    for s in sess_data_files:
        if re.search(sess_key['session_id'], s):
            return s
    return None


def split_list(arr, chunk_size):
    for s in range(0, len(arr), chunk_size):
        yield arr[s:s+chunk_size]
