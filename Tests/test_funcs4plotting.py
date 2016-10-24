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

@slow
def test_save_plot(test_data):
    map_plot(test_data[:,:,0])

    save_plot()
    filename_0 = "myplot.png"

    save_plot(title="test_1")
    filename_1 = "test_1.png"

    save_plot(title="test_2", location=wd)
    filename_2 = os.path.join(wd, "test_2.png")

    save_plot(title="test_3", location="new_folder")
    filename_3 = os.path.join("new_folder", "test_3.png")

    save_plot(title="test_4", location=wd, extensions=["pdf", "png"])
    filename_4 = os.path.join(wd, "test_4.png")
    filename_5 = os.path.join(wd, "test_4.pdf")

    # Test plot has been created, and then remove it"
    filenames = [filename_0, filename_1, filename_2, 
                filename_3, filename_4, filename_5]
    for filename in filenames:
        assert os.path.isfile( filename ), "Failed to create {file}".format(file=filename)
        os.remove(filename)

    os.rmdir("new_folder")

    return
    





