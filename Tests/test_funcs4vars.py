from ..bpch2netCDF import *
import logging
import pytest
logging.basicConfig(filename='test.log',level=logging.DEBUG)                    
logging.info('Starting funcs4GEOSC test.') 

noCTM = pytest.mark.skipif(                                                     
    not pytest.config.getoption("--remake_ctm"),                                
    reason="need --runslow option to run"                                       
)   

@noCTM
def test_convert_to_netCDF():
    return

def test_get_folder():
    return





logging.info('funcs4GEOSC test complete')
