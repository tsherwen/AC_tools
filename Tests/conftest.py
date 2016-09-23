#####
# This file contains the settings used for py.test
#####

import pytest
import logging

#aFORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"            
FORMAT = "%(filename)s:%(lineno)s - %(funcName)s() : %(message)s"

logging.basicConfig(filename='test.log',level=logging.DEBUG, format=FORMAT)
                                             


def pytest_addoption(parser):
    parser.addoption("--remake_ctm", action="store_true",
        help="remake ctm.nc tests")
