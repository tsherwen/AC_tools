from ..funcs4plotting import *
import logging
import pytest
logging.basicConfig(filename='test.log',level=logging.DEBUG)

slow = pytest.mark.skipif(                                                      
    not pytest.config.getoption("--slow"),                                      
    reason="need --slow option to run"                                          
)                                                                               
        
wd='test_files'
@pytest.fixture()
def test_data():
    from ..funcs4GEOSC import get_GC_output 
    test_data = get_GC_output('test_files', species='O3') 
    return test_data

@slow
def test_map_plot(test_data):
    logging.info("begining test")
    map_plot( test_data[:,:,0] )
    map_plot( test_data[:,:,0].T)
    with pytest.raises(AssertionError):
        map_plot( test_data[0,:,:], wd=wd )
        map_plot( None )
    return





