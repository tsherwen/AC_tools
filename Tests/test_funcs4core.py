from ..funcs4core import *
import logging
import pytest
logging.basicConfig(filename='test.log',level=logging.DEBUG)                    

noCTM = pytest.mark.skipif(                                                     
    not pytest.config.getoption("--remake_ctm"),                                
    reason="need --runslow option to run"                                       
)   

@noCTM
def test_convert_to_netCDF():
    return

def test_gchemgrid():
    # Test fail for no inputs
    with pytest.raises(Exception):
        gchemgrid()

    # Test dictionary return
    arr = gchemgrid('c_km_geos5_r')
    assert isinstance( arr, np.ndarray), 'item is not a numpy array'
    
    dic = gchemgrid( rtn_dict=True )
    assert isinstance( dic, dict), 'Dictionary return failed.'
    return
    




