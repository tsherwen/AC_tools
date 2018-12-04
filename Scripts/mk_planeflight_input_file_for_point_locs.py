#!/usr/bin/python
# --- Packages
import numpy as np
from time import gmtime, strftime
import datetime as datetime
import time
import glob
import pandas as pd
import sys
from . import AC_tools as AC

def main(filename=None, LAT_var='LAT', LON_var='LON', \
        PRESS_var='PRESS', loc_var='TYPE', Username='Tomas Sherwen', \
        slist=None, Extra_spacings=False, freq='H', start_year=2014, \
        end_year=2016, GC_ver_above_v12=True, GC_ver='v12.0.0', \
        debug=False):
    """
    Mk planeflight input files for GEOS-Chem ND40 diagnostic from a csv file 
    contains sites of interest (TYPE) and their descriptors (LON, LAT, PRESS)
    
    Parameters
    -------
    freq (str): output frequency of diagnostic ('H' or 'D')
    start_year,end_year (int): start and end year to output.
    wd (str): the working (code) directory to search for files in
    loc_var (str): name for (e.g. plane name), could be more than one. 
    LAT_var, LON_var, PRESS_var (str): name for pressure(HPa),lat and lon in df
    Username (str): name of the programme's user
    Extra_spacings (boolean): add extra spacing? (needed for large amounts of 
        output, like nested grids)
    Notes
    -----
     - command line call with the following:
    "python mk_planeflight_input_file_for_point_locs.py <name of dat file>"
    ( default file provided (Planeflight_point_sites.csv) to edit)
     - This programme can be used to produce files to output data for 
    observational sites (e.g. goverment air quality sites)
     - This programme makes the planeflight*.dat files required to output 
    for specific locations and times in the model. 
    - The deafult settting is for hourly output.
    """
    print('filename:{}'.format(filename)) 
    # ---  Local settings 
    # file of locations? (from 1st argument of command line)
    if isinstance(filename, type(None)):
        if __name__ == "__main__":
            filename = sys.argv[1]
        else:
            filename = 'Planeflight_point_sites.csv'    

    # set species list to output if not provided
    if isinstance(slist, type(None)):
        # Which (halogen) code version is being used?
        #ver = '1.6' # Iodine simulation in v9-2
        #ver = '2.0' # Iodine + Bromine simulation
        ver = '3.0' # Cl-Br-I simulation 
        # Get Variables to output (e.g. tracers, species, met values )
        slist = AC.pf_var('slist', ver=ver, fill_var_with_zeroes=True )

    # --- Read in site Detail
    # ( must contain LAT, LON, PRESS, TYPE (name of sites) )
    
    LOCS_df = pd.read_csv( filename )
    vars_ = ['LAT', 'LON', 'PRESS', 'TYPE' ]
    LAT, LON, PRESS, TYPE = [ LOCS_df[i].values for i in vars_ ]

    # --- Make DataFrame of locations 
    dates = pd.date_range( datetime.datetime(start_year, 1, 1), \
        datetime.datetime(end_year, 12, 31, 23), freq='H' )
    # for each location make a DataFrame, then conbime
    dfs =[]
    for n, type_ in enumerate( TYPE ):
        # dictionary of data
        nvar = len(dates)
        d = {
        'datetime':dates,'LAT':[LAT[n]]*nvar, 'LON':[LON[n]]*nvar,
        'TYPE':[TYPE[n]]*nvar, 'PRESS':[PRESS[n]]*nvar}
        dfs += [ pd.DataFrame(d, index=np.arange(nvar)+(n*1E6)) ]
    # combine all TYPE (sites) and sort by date
    df = pd.concat(dfs).sort_values('datetime',ascending=True)
    
    # --- Print out files 
    if (GC_ver == 'v12.0.0') or GC_ver_above_v12:
        AC.prt_PlaneFlight_files_v12_plus(df=df, slist=slist, Extra_spacings=Extra_spacings)
    else:
        AC.prt_PlaneFlight_files(df=df, slist=slist, Extra_spacings=Extra_spacings)


if __name__ == "__main__":
    main()

