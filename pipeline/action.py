'''
Schema of session information.
'''
import datajoint as dj
from pipeline import reference, subject

schema = dj.schema(dj.config.get('database.prefix', '') + 'hg2015_action')


@schema
class Weighing(dj.Manual):
    definition = """
    -> subject.Subject
    ---
    weight_before: float   # in grams
    weight_after: float    # in grams
    """


@schema
class SubjectWhiskerConfig(dj.Manual):
    definition = """
    -> subject.Subject
    ---
    -> reference.WhiskerConfig
    """
