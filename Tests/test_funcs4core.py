from ..funcs4core import *
import logging
import pytest
logging.basicConfig(filename='test.log',level=logging.DEBUG)                    

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
    




