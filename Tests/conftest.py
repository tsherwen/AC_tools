#####
# This file contains the settings used for py.test
#####

import pytest
import logging
import os
# temporality restore urllib2 as urllib not installed
#import urllib.request, urllib.error, urllib.parse
import urllib2

#FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"            
FORMAT = "%(filename)s:%(lineno)s - %(funcName)s() : %(message)s"
test_file_dir = '../data'      

logging.basicConfig(filename='test.log',level=logging.DEBUG, format=FORMAT)
logging.getLogger().setLevel(logging.DEBUG)                                     

                                             


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true",
        help="remake ctm.nc tests")

def pytest_configure():

    # Make sure we are in the correct folder for the test.
    dirname = os.path.split(os.getcwd())[1]
    if not dirname=='Tests':
        pytest.exit("Not running in the Tests folder!")


    # Make sure the data is downloaded
    from ..Scripts import get_data_files
                                                                                
    return           
