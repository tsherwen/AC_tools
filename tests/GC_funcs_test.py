#import netCDF4
from .. import io
from .. import GC_funcs
import logging
logging.basicConfig(filename='test.log', level=logging.DEBUG)
logging.info('Started GC_funcs test')

#pytest_plugins = ['pytest_profiling']

 
#GC_funcs.get_variables( netCDF_file )

netCDF_file = io.open_netCDF(filename='test_files/ctm.nc')                             

def test_get_surface_area():
    logging.info("Begining test.")
    surface_area = GC_funcs.get_surface_area( netCDF_file )
    assert len(surface_area.shape) == 2
    logging.info("Test complete.")
    return

def test_get_variables():
    logging.info("Begining test.")
    variables = GC_funcs.get_variables( netCDF_file, show=False )
    assert isinstance(variables, list) == True
    assert variables[0] == "ACETSRCE__ACETbg"
    assert variables[-1] == "WETDLS_S__SO4s"
    logging.info("Test complete.")
    return

def test_get_variable_data():
    variable = 'IJ_AVG_S__O3'
    new_data = GC_funcs.get_variable_data( netCDF_file, variable )
    assert len(new_data.shape) == 3 

    variable2 = 'fail'
    try:
        new_data2 = GC_funcs.get_variable_data( netCDF_file, variable2 )
        test = False
    except IOError:
        test = True
    assert test, 'Failed catching a fake variable'
    return


# Currently dont have land map in the bpch or netCDF so cannot check this is working.
#def test_get_land_map():
#    data = GC_funcs.get_land_map(netCDF_file)
#    assert( len(land_map.shape) == 3 ), 'land map shape is wrong.'
#    return

def test_get_air_mass():
    air_mass = GC_funcs.get_air_mass(netCDF_file)
    assert( len(air_mass.shape) == 3 ), 'air mass shape is wrong.'
    return

def test_get_trop_time():
    trop_time = GC_funcs.get_trop_time(netCDF_file)
    assert(len(trop_time.shape) > 0), 'No shape for the trop time.'
#    assert( len(air_mass.shape) == 4 ), 'trop time shape is wrong.'
    return

def test_get_species_rmm():
    species_rmm = GC_funcs.get_species_rmm( 'O3' )
    assert species_rmm == 47.98474386, \
        "Species rmm is not working for defined species (O3)"
    species_rmm = GC_funcs.get_species_rmm( 'bob' )
    assert species_rmm == 1.0, \
        "Species rmm is not working for undefiend species"
    return

def test_get_tropospheric_burden():
    import numpy as np
    tropospheric_burden = GC_funcs.get_tropospheric_burden( netCDF_file, 'O3' )
    assert isinstance( tropospheric_burden, float), 'The tropospheric burden is not a float.'
    return

def test_get_species_data():
    species_data = GC_funcs.get_species_data( netCDF_file, 'O3' )
    assert( len(species_data.shape) == 3), 'species data shape is wrong'
    return

def test_get_volume():
    volume = GC_funcs.get_volume( netCDF_file )
    assert( len(volume.shape) == 3), 'Volume data is wrong.'
    return

def test_get_tropospheric_PL():
    O3_production = GC_funcs.get_tropospheric_PL( netCDF_file, "PO1", 48.0 )
    assert( len(O3_production.shape) == 3), 'tropsopheric pl is wrong.'
    return

def test_get_tropospheric_total_pl():
    O3_production = GC_funcs.get_tropospheric_total_PL( netCDF_file, "PO1", 48.0 )
    assert isinstance(O3_production, float), "total tropsopheric pl is not a float."
    return
 
def test_get_drydep():
    O3_drydep = GC_funcs.get_drydep( netCDF_file, "O3")
    assert (len(O3_drydep.shape) == 2), 'drydep shape is wrong.'
    return

def test_get_annual_drydep():
    O3_drydep = GC_funcs.get_annual_drydep( netCDF_file, "O3")
    assert isinstance(O3_drydep, float), 'Total annual drydep is not a float.'


logging.info('Finished GC_funcs test')
