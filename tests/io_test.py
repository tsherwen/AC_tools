"""
Unit test to test all functions in AC_tools.io
Expected to be slow due to IO and netCDF creation.
"""


from .. import io
import logging

# Make sure the files are downloaded

logging.basicConfig(filename='test.log', level=logging.DEBUG)

logging.info('Started io test')

def test_open_netCDF():
    """
    Create a netCDF file from bpch and confirm that all variables inside are equivelent.
    """
       
    import netCDF4
    import os

    assert os.path.isfile("test_files/ctm.nc"), "Couldnt find the test ctm.nc file.  \
    have you run the setup.py ?"

    test_netCDF_data = netCDF4.Dataset(filename='test_files/ctm.nc')
    test_bpch_data = io.open_netCDF(folder='test_files',filename='temp.nc')
    for item1, item2 in zip(test_netCDF_data.variables, test_bpch_data.variables):
        assert item1 == item2

    # Clean up
    os.remove('test_files/temp.nc')
    return

logging.info('Finished io test')
