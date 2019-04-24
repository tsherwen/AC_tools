from ..core import *
import logging
import pytest
import numpy as np
logging.basicConfig(filename='test.log', level=logging.DEBUG)

data_dir = '../data'


def test_get_sigfig():
    print(get_sigfig(100, 1))
    assert(get_sigfig(3.22294, 1) == 3)
    assert(get_sigfig(3.29294, 2) == 3.3)
    assert(get_sigfig(3.29294, 3) == 3.29)
    return


def test_get_sigfig_big():
    #    assert(get_sigfig(3.22294E20, 0)==3E20)
    assert(get_sigfig(3.22294E20, 1) == 3E20)
    assert(get_sigfig(3.29294E20, 2) == 3.3E20)
    assert(get_sigfig(3.29294E20, 3) == 3.29E20)


def test_get_sigfig_small():
    assert(get_sigfig(3.22294E-20, 1) == 3E-20)
    assert(get_sigfig(3.29294E-20, 2) == 3.3E-20)
    assert(get_sigfig(3.29294E-20, 3) == 3.29E-20)
    return

# def test_get_sigfig_negative():
#    assert(get_sigfig(-3.22294, 1)==-3.2)
#    assert(get_sigfig(-3.29294, 2)==-3.29)
#    assert(get_sigfig(-3.29294, 3)==-3.293)
#    return


def test_get_scientific_number():
    assert(get_scientific_number(3.22294, 1, string=True) == "3.2")
    assert(get_scientific_number(3.29294E10, 2, string=True) == "3.29E10")
    assert(get_scientific_number(-3.29294E-10, 3, string=True) == "-3.293E-10")
    return


def test_gchemgrid():
    # Test fail for no inputs
    with pytest.raises(Exception):
        gchemgrid()

    # Test dictionary return
    arr = gchemgrid('c_km_geos5_r')
    assert isinstance(arr, np.ndarray), 'item is not a numpy array'

    dic = gchemgrid(rtn_dict=True)
    assert isinstance(dic, dict), 'Dictionary return failed.'
    return


def test_get_dims4res_4x5():
    # test the dictionary returns
    res = get_dims4res('4x5')
    assert (res == (72, 46, 47)), '4x5 res lookup failed'


def test_get_dims4res_2x25():
    # test the dictionary returns
    res = get_dims4res('2x2.5')
    assert (res == (144, 91, 47)), '2x2.5 res lookup failed'


def test_get_dims4res_2D():
    res_2D = get_dims4res('4x5', just2D=True)
    assert (res_2D == (72, 46)), "2D lookup failed"


def test_get_dims4res_dict():
    test_dict_2d = get_dims4res(r_dims=True)
    assert (test_dict_2d[72, 46, 47] == '4x5'), 'dict lookup failed'


def test_get_latlonalt4res_default():

    # Test default
    (lon, lat, alt) = get_latlonalt4res()
    print("lon = ", lon)
    assert len(lat) == 46, 'The default latitude is wrong'
    assert len(lon) == 72, 'The default longitude is wrong'
    assert len(alt) == 47, 'The default altidure is wrong'
