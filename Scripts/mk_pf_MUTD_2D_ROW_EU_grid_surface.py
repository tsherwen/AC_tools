#!/usr/bin/python
# --- Packages
import numpy as np
from time import gmtime, strftime
import time
import glob
import AC_tools as AC
import sys
import pandas as pd


def main(filename=None, LAT_var='LAT', LON_var='LON',
         PRESS_var='PRESS', loc_var='TYPE', Username='Tomas Sherwen',
         slist=None, Extra_spacings=True,
         freq='H', start_year=2015, end_year=2017, debug=False):
    """
    This programme makes the planeflight*.dat files required to output for 
    an entire grid (e.g. Europe.)

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
    res = '0.25x0.3125'  # res='0.5x0.666'
    tpwd = AC.get_dir('tpwd')
    #start_year, end_year = 2006,2007
    # location of files? (from 1st argument of command line)
    if isinstance(filename, type(None)):
        filename = sys.argv[1]
    else:
        filename = 'EU_GRID_{}.dat'.format(res)

    # set species list to output if not provided
    if isinstance(slist, type(None)):
        # Which (halogen) code version is being used?
        # ver = '1.6' # Iodine simulation in v9-2
        # ver = '2.0' # Iodine + Bromine simulation
        ver = '3.0'  # Cl-Br-I simulation
        # Get Variables to output (e.g. tracers, species, met values )
        slist = AC.pf_var('slist', ver=ver, fill_var_with_zeroes=True)

    # --- Read in site Detail
    location = AC.readin_gaw_sites(filename, all=True)
    numbers, locs, lats, lons, pres = [location[:, i] for i in range(5)]
    lats, lons, pres = [np.float64(i) for i in (lats, lons, pres)]
    locs = np.array(locs)
    print((lats[0:4]))

    # --- Set Variables
    # slist = pf_var( 'slist_v9_2_NREA_red_NOy', ver=ver )#'slist_v9_2_NREA_red' )
    #slist = pf_var( 'slist_ClearFlo', ver=ver )
    slist = AC.pf_var('slist_PEN_WEY_LEI', ver=ver)
    # kludge additional of ethanol
    slist = slist + ['TRA_86']
    nvar = len(slist)
    # dates
    dates = pd.date_range(datetime.datetime(start_year, 1, 1),
                          datetime.datetime(end_year, 12, 31, 23), freq='H')
    # dictionary of variables
    d = {
        'datetime': dates, 'LAT': [LAT[n]]*nvar, 'LON': [LON[n]]*nvar,
        'TYPE': [TYPE[n]]*nvar, 'PRESS': [PRESS[n]]*nvar}
    df = pd.DataFrame(d)

    # --- Print out files
    AC.prt_PlaneFlight_files(df=df, slist=slist, Extra_spacings=Extra_spacings)


if __name__ == "__main__":
    main()
