from ..funcs4GEOSC import *
import logging
import pytest
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(filename='test.log',level=logging.DEBUG, format=FORMAT)
logging.info('Starting funcs4GEOSC test.')

wd = 'Test_files/GC_run'

# Option to remake CTM files - Not recomended for normal tests - slow.
noCTM = pytest.mark.skipif(
    not pytest.config.getoption("--remake_ctm"),
    reason="need --runslow option to run"
)

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
