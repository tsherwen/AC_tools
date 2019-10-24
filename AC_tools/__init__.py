# compatibility with both python 2 and 3
from __future__ import print_function
import numpy as np
import sys
# AC_tools modules
from . AC_time import *
from . core import *
from . generic import *
from . GEOSChem_bpch import *
from . GEOSChem_nc import *
from . KPP import *
from . mask import *
from . planeflight import *
from . plotting import *
from . SMVGEAR import *
from . variables import *
# include the redundent files for now
from . obsolete.plotting_REDUNDANT import *
from . obsolete.variables_REDUNDANT import *

"""
AC_tools is a module of functions started by Tomas, and contributed to by others in the Evans' group, and hopefully maintained by the Group.
To access the help, from python or ipython, type help(AC_tools) to get general help
To get more detailed help from a module for example, type help(AC_tools.AC_time.py)
If you find missing documentation any thing is unclear in any of this, please request a git push to github.
"""

# Setup logging for module
import logging
level = logging.DEBUG
FORMAT = "%(levelname)8s - %(message)s   @---> %(filename)s:%(lineno)s  %(funcName)s()"
logging.basicConfig(filename='AC_tools.log', filemode='w', level=level,
                    format=FORMAT)
logging.getLogger().setLevel(level)

# Import submodules here for easier access
