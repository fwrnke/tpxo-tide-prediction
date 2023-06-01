import os  # noqa
from .predict_tide import (  # noqa
    read_parameter_file,
    tide_predict,
    write_tides
    )

__version__ = open(os.path.join(os.path.dirname(__file__), '../VERSION')).read().rstrip('\n')
