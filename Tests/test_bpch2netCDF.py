from ..bpch2netCDF import *
import logging
import pytest
import os

noCTM = pytest.mark.skipif(                                                     
    not pytest.config.getoption("--remake_ctm"),                                
    reason="need --runslow option to run"                                       
)   

wd='Test_files/GC_run'

@noCTM
def test_convert_to_netCDF():
    logging.info("beginning test")
    import filecmp
    # Recreate a ctm.nc file and confirm it is the same
    logging.debug("Creating the temp netCDF file")
    convert_to_netCDF(folder=wd, filename='temp.nc')
    datafile = os.path.join(wd, 'ctm.nc')
    tempfile = os.path.join(wd, 'temp.nc')

    logging.debug("Comparing the temp netCDF file to the origional")
    assert filecmp.cmp(datafile, tempfile), \
        'bpch converter failed to replicate the origional file.'

    os.remove(tempfile)
    logging.info("test complete")
    return

def test_get_folder():
    logging.info("beginning test")
    folder = get_folder(wd)
    assert isinstance(folder, str), "The folder is not a string"
    import os.path
    assert os.path.exists(folder), "Cannot find the test folder"
    logging.info("test complete")
    return





