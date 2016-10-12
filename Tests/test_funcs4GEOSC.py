from ..funcs4GEOSC import *
import logging
import pytest
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(filename='test.log',level=logging.DEBUG, format=FORMAT)
logging.info('Starting funcs4GEOSC test.')

wd = 'test_files'


#1.04
def test_get_surface_area():
    arr = get_surface_area()
    assert isinstance( arr, np.ndarray), 'Surface area not a numpy array.'
    assert (len(arr.shape)==3), 'Surface area has too many dimensions.'
    return

#1.05
def get_land_map():
    arr = get_land_map()
    assert isinstance( arr, np.ndarray), 'Land map not a numpy array.'
    assert (len(arr.shape) == 3), 'Land map has too many dimensions.'
    return



#1.07
def test_get_air_mass_np():
    logging.info("Begining test.")
#    from funcs4GEOSC import get_air_mass_np
    arr = get_air_mass_np(wd=wd)
    assert isinstance( arr, np.ndarray), 'Air mass array is not a numpy array'
    logging.info("Test complete.")
    return

#1.22
def test_get_GC_output():
    logging.info("Beginning test.")
    arr = get_GC_output(wd=wd, species='O3', category='IJ_AVG_S')
    assert isinstance( arr, np.ndarray), 'GC output is not a numpy array'
    logging.info("Test complete.")
    return

logging.info('funcs4GEOSC test complete.')   
