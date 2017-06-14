#!/usr/bin/python
"""
Check GEOS-Chem met fields for zero in arrays.

NOTES:
 - credit: mje, adapted by tms

"""
# Import modules
from netCDF4 import Dataset
import numpy as np
import glob as glob
import matplotlib.pyplot as plt

# ---- Setup 
#
# For what years and months? ("*" for all )
years='*'#'2014'
months='*'
# on what grid?
grid = 'GEOS_0.25x0.3125_ch' # GEOS_0.25x0.3125_eu
# where is the data?
data_root = '/work/data/GEOS'
# meterology
met='GEOS_FP'
# add extension?
ext = '*.I3*' # '*.nc*'  # I3 files or all?
# directory to check?
dir='/'.join( (data_root, grid, met, years, months, ext) )
# get files
files=glob.glob(dir)
files.sort()
# print detail to screen?
prt_detail=False

# --- loop files and test for zero values
print(files)
counter=0
for file in files:
	# Debug? - print file accessed... 
#	print file
	# Open NetCDF as dataset ("d" )	
	with Dataset(file) as d:
	# loop keys
		for key in list(d.variables.keys()):
			# check on a per field basis
			field=d[key][:]
			# If multi-dimensional
			if (len(field.shape) > 1):
				# print zero fields
				if (field.min() == 0.):
					# print to screen
					print(file, key, field.shape)
					# Print detail?
					if prt_detail:
						for j in np.arange(0,8):
							for k in np.arange(0,72):
								if (field[j,k,:,:].min() == 0.):
									print(key,j,k,file)
