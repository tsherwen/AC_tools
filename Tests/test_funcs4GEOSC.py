from ..funcs4GEOSC import *
import logging
logging.basicConfig(filename='test.log',level=logging.DEBUG)
logging.info('Starting funcs4GEOSC test.')

wd = 'test_files/GC_run'

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
