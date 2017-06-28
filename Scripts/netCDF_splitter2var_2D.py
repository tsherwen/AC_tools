#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Split off 2D variable from file with other variables

Notes
---- 
 - based on software carpentary example. 
http://damienirving.github.io/capstone-oceanography/03-data-provenance.html
"""
# Modules to import
from netCDF4 import Dataset
import numpy as np
import pylab as pl
import calendar
# add extra's for copied function... 
import os, sys, argparse
import datetime

# --- verbose and debug settings for script main call 
VERBOSE=False
DEBUG=False


def main( filename=None, VarName='OLSON', verbose=False, debug=False ):
    """
    Driver to split off variables
    """
    # Get the file name and location
    wd, fn = get_file_loc_and_name()
    # name output file if name not given
    if isinstance( filename, type(None) ):
        filename = wd.split('/')[-2]
    if debug:
        print((wd, fn, filename))
    inFile = wd+'/'+fn

    # Set output name
    outfile_name = inFile+'.out'
    
    # Read input data 
    VarData, input_DATA = read_data(inFile, VarName=VarName)
    
    # Set values?
#    print type(VarData)
#    print [ (i.shape, i.mean(), i.min(), i.max()) for i in VarData]
#    VarData[VarData>1] = 1
#    print [ (i.shape, i.mean(), i.min(), i.max()) for i in VarData]
    
    # --- Write the output file
    outfile = Dataset(outfile_name, 'w', format='NETCDF4')
    set_global_atts(input_DATA, outfile)
    copy_dimensions(input_DATA, outfile)
    copy_variables(input_DATA, outfile, VarName=VarName)
    # overwite data
    outfile[VarName][:] = VarData
    # Close file
    outfile.close()


def get_file_loc_and_name( ):
    """ Get file location and name """

    # Use command line grab function
    import sys
    
    # Get arguments from command line
    wd = sys.argv[1]
    fn = sys.argv[2]

    return wd, fn 
    

def copy_dimensions(infile, outfile):
    """ 
    Copy the dimensions of the infile to the outfile
    """    
    for dimName, dimData in iter(list(infile.dimensions.items())):
        outfile.createDimension(dimName, len(dimData))


def copy_variables(infile, outfile, VarName='OLSON'):
    """
    Create variables corresponding to the file dimensions 
    by copying from infile
    """
    # Get vars
    var_list = ['lon', 'lat', 'time']
    # Also consider LANDMAP value
    var_list+=[VarName]
    # Now loop        
    for var_name in var_list:
        varin = infile.variables[var_name]
        outVar = outfile.createVariable(var_name, varin.datatype, 
                                        varin.dimensions, 
                                        )
        outVar[:] = varin[:]
            
        var_atts = {}
        for att in varin.ncattrs():
            if not att == '_FillValue':
                var_atts[att] = eval('varin.'+att) 
        outVar.setncatts(var_atts)


def read_data(ifile, VarName='OLSON'):
    """ 
    Read data from ifile corresponding to the VarName 
    """
    input_DATA = Dataset(ifile)
    VarData = input_DATA.variables[VarName][:]

    return VarData, input_DATA


def set_global_atts(infile, outfile):
    """Set the global attributes for outfile.
        
    Note that the global attributes are simply copied from infile.
    """
        
    global_atts = {}
    for att in infile.ncattrs():
        global_atts[att] = eval('infile.'+att)  
        
    # set attributes
    outfile.setncatts(global_atts)


if __name__ == "__main__":
    main( verbose=VERBOSE, debug=DEBUG )
