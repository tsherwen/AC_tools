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

def test_get_dims4res():
    # test the dictionary returns
    test_dict = get_dims4res( r_dims=True )
    assert (test_dict[72,46,47] == '4x5'), '4x5 res lookup failed'

    test_dict_2d = get_dims4res( r_dims=True, just2D=True)
    assert (test_dict_2d[72,46] == '4x5'), '4x5 res lookup failed'
    




