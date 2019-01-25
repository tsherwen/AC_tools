from ..bpch2netCDF import *
import logging
import pytest
import os
#import urllib.parse, urllib.error
# temporality restore urllib2 as urllib not installed
from urllib import urlopen
import urllib2

slow = pytest.mark.skipif(
    not pytest.config.getoption("--slow"),
    reason="need --slow option to run"
)

test_file_dir = '../data'

# def setup_function(function):
#    """
#    Downloads all the test dataset files using rsync.
#    """
#
#    test_files = ['test.nc', 'test.bpch','tracerinfo.dat','diaginfo.dat']
#
#    url_base =  'http://atmosviz1.york.ac.uk/~bn506/data/AC_tools/'
#    test_file_dir = 'test_files'
#
#
#    if not os.path.exists(test_file_dir):
#        os.makedirs(test_file_dir)
#
#    for file_name in test_files:
#        file_path = os.path.join(test_file_dir, file_name)
#        if not os.path.isfile(file_path):
#            my_file = open(file_path, 'wb')
#            logging.debug(file_name + " not found. Downloading now.")
#            url = url_base + file_name
#            file_data = urllib2.urlopen( url ).read()
#            my_file.write(file_data)
#            my_file.close()
#
#            logging.debug(file_name + " downloaded.")
#
#    return


def file_comparison(file_1, file_2):
    file_1_data = open(file_1, 'r')
    file_2_data = open(file_2, 'r')
    if file_1_data.read() == file_2_data.read():
        same = True
    else:
        same = False
    return same


@slow
def test_convert_to_netCDF():
    logging.info("beginning test")
    # Recreate a ctm.nc file and confirm it is the same
    logging.debug("Creating the temp netCDF file")
    convert_to_netCDF(folder=test_file_dir, bpch_file_list=[
                      'test.bpch'], remake=True, filename='test.nc')
    datafile = os.path.join(test_file_dir, 'ctm.nc')
    testfile = os.path.join(test_file_dir, 'test.nc')

    logging.debug("Comparing the temp netCDF file to the origional")
    assert file_comparison(datafile, testfile), \
        'bpch converter failed to replicate the original file.'

    os.remove(testfile)
    logging.info("test complete")
    return


def test_get_folder():
    logging.info("beginning test")
    folder = get_folder(test_file_dir)
    assert isinstance(folder, str), "The folder is not a string"
    assert os.path.exists(folder), "Cannot find the test folder"
    logging.info("test complete")
    return
