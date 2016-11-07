from ..funcs4core import *
import logging
import pytest
import numpy as np
logging.basicConfig(filename='test.log',level=logging.DEBUG)                    

data_dir = '../data'

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


def test_get_latlonalt4res_default():

    # Test default
    lon, lat, alt = get_latlonalt4res()
    print type(lon)
    assert isinstance(lon, np.ndarray), 'The lon is not correct.'
    assert isinstance(lat, np.ndarray), 'The lat is not correct.'
    assert isinstance(alt, np.ndarray), 'The alt is not correct.'
    
def test_get_latlonalt4res_res():
    lon, lat, alt = get_latlonalt4res(res='4x5')
    assert isinstance(lon, np.ndarray), 'The lon is not correct.'
    assert isinstance(lat, np.ndarray), 'The lat is not correct.'
    assert isinstance(alt, np.ndarray), 'The alt is not correct.'
    # Test give_resolution

def test_get_latlonalt4res_wd():
    lon, lat, alt = get_latlonalt4res(wd=data_dir)
    assert isinstance(lon, np.ndarray), 'The lon is not correct.'
    assert isinstance(lat, np.ndarray), 'The lat is not correct.'
    assert isinstance(alt, np.ndarray), 'The alt is not correct.'
    # Test give wd



