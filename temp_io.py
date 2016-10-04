

from GC_tools.io import open_netCDF
import netCDF4
netCDF_data = open_netCDF('tests/test_files')

ben_data = netCDF_data.variables['IJ_AVG_S__BEN']
netCDF_data.createDimension('time')

print ben_data
