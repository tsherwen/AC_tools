from ..bpch2netCDF import *
import logging
import pytest
import os
import wget
import filecmp

slow = pytest.mark.skipif(                                                     
    not pytest.config.getoption("--slow"),                                
    reason="need --slow option to run"                                       
)   

test_file_dir = 'test_files'

def setup_function(function):
    """                                                                                    
    Downloads all the test dataset files using rsync.                                      
    """

    test_files = ['ctm.nc', 'test.bpch','tracerinfo.dat','diaginfo.dat']

    url_base =  'http://atmosviz1.york.ac.uk/~bn506/data/GC_funcs_test/test_files/'
    test_file_dir = 'test_files'


    if not os.path.exists(test_file_dir):
        os.makedirs(test_file_dir)

    for file_name in test_files:
        file_path = os.path.join(test_file_dir, file_name)
        if not os.path.isfile(file_path):
            logging.debug(file_name + " not found. Downloading now.")
            url = url_base + file_name
            filename = wget.download(url, out=test_file_dir)
            logging.debug(file_name + " downloaded.")
        
            
                                                                     
#    if not os.path.isfile("test_files/ctm.nc"):
#        logging.debug("Downloading test files")
#        os.
#        os.system("""(cd test_files && wget -r -nH -nd -np -R --no-parent --reject "index.html*" http://atmosviz1.york.ac.uk/~bn506/data/GC_funcs_test/test_files/)""")
#        logging.debug("Download of test files complete")
#        # Consider making the above more pythonic.
    return


@slow
def test_convert_to_netCDF():
    logging.info("beginning test")
    # Recreate a ctm.nc file and confirm it is the same
    logging.debug("Creating the temp netCDF file")
    convert_to_netCDF(folder=test_file_dir, bpch_file_list=['test.bpch'], remake=True)
    datafile = os.path.join(test_file_dir, 'ctm.nc')
    testfile = os.path.join(test_file_dir, 'test.nc')

    logging.debug("Comparing the temp netCDF file to the origional")
    assert filecmp.cmp(datafile, testfile), \
        'bpch converter failed to replicate the origional file.'

    os.remove(datafile)
    logging.info("test complete")
    return

def test_get_folder():
    logging.info("beginning test")
    folder = get_folder(test_file_dir)
    assert isinstance(folder, str), "The folder is not a string"
    assert os.path.exists(folder), "Cannot find the test folder"
    logging.info("test complete")
    return





