"""
AC_tools is a module of functions started by Tomas, and contributed to by others in the Evans' group, and hopefully maintained by the Group.
To access the help, from python or ipython, type help(AC_tools) to get general help
To get more detailed help from a module for example, type help(AC_tools.funcs4time.py)
If you find missing documentation any thing is unclear in any of this, please request a git push to github.
"""

import logging
FORMAT = "%(levelname)8s - %(message)s   @---> %(filename)s:%(lineno)s  %(funcName)s()"
logging.basicConfig(filename='AC_tools.log', filemode='w',level=logging.DEBUG,
                format=FORMAT)
logging.getLogger().setLevel(logging.DEBUG)

import numpy as np
# Import submodules here for easier access
from . funcs4GEOSC import *
from . funcs4core import *
from . funcs4generic import *
from . funcs4pf import *
from . funcs4time import *
from . funcs_vars import *
from . funcs4plotting import *
