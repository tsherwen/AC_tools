#!/usr/bin/python
# --- Packages
import numpy as np
from time import gmtime, strftime
import time
import glob
from . import AC_tools as AC
import sys
import pandas as pd


def main(filename=None, slist=None, time_str='%H:%M', date_str='%d%m%Y',
         LAT_var='LAT', LON_var='LON', PRESS_var='PRESS', loc_var='TYPE',
         Username='Tomas Sherwen', GC_ver='v12.0.0', GC_ver_above_v12=True):
    """
    This programme makes the planeflight*.dat files required to output for 
    specific locations and times in the model (e.g. along a ship cruise or 
    airborne campaign track).

    Parameters
    -------
    time_str (str): format of time columne in input .csv file

    Returns
    -------

    Notes
    -----
     - This programme can be used to produce files to output data for ship 
    and aricraft campaigns
    """
    # ---  Local settings
    # location of files? (from 1st argument of command line)
    if isinstance(filename, type(None)):
        filename = sys.argv[1]
    # set species list to output if not provided
    if isinstance(slist, type(None)):
        # Which (halogen) code version is being used?
        # ver = '1.6' # Iodine simulation in v9-2
        # ver = '2.0' # Iodine + Bromine simulation
        ver = '3.0'  # Cl-Br-I simulation
        # Get Variables to output (e.g. tracers, species, met values )
        slist = AC.pf_var('slist', ver=ver, fill_var_with_zeroes=True)

    # --- Read in a .csv file of aircraft/cruise data etc...
    # (must contrain DATE, TIME, LAT, LON, PRESS, (hPa), TYPE (code ) )
    # (the tag for output (referred to as TYPE in GC) can be upto 4 chars )
    LOCS_df = pd.read_csv(filename)
    # or multiple input files
#    mk_Campaign_csv_from_multiple_csv_files()
#    vars_ = ['LAT', 'LON', 'PRESS', 'TYPE' ]
#    LAT, LON, PRESS, TYPE = [ LOCS_df[i].values for i in vars_ ]
#    nvar = len( slist )
    df = LOCS_df

    # --- Print out files
    if (GC_ver == 'v12.0.0') or GC_ver_above_v12:
        AC.prt_PlaneFlight_files_v12_plus(df=df, slist=slist, Extra_spacings=Extra_spacings,
                                          LON_var=LON_var, LAT_var=LAT_var, PRESS_var=PRESS_var, loc_var=loc_var,
                                          Username=Username)
    else:
        AC.prt_PlaneFlight_files(df=df, slist=slist, Extra_spacings=Extra_spacings,
                                 LON_var=LON_var, LAT_var=LAT_var, PRESS_var=PRESS_var, loc_var=loc_var,
                                 Username=Username)


if __name__ == "__main__":
    main()
