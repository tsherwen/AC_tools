#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Functions for use with the Harvard-NASA Emissions Component (HEMCO)

Use help(<name of function>) to get details on a particular function.

"""
# - Required modules:
# I/O / Low level
import os
import sys
import glob
import pandas as pd
import logging
# Math/Analysis
import numpy as np
# Time
import time
import datetime as datetime


def rm_eruptive_volcancos_from_files(dates=None, sdate=None, edate=None,
                                     folder=None, output_folder=None):
    """
    Remove the eruptive volcanos from the input HEMCO emission files
    """
    # Set dates to update files as NASA ATom4, unless others provided
    if isinstance(dates, type(None)):
        if isinstance(sdate, type(datetime.datetime)):
            sdate = datetime.datetime(2018, 1, 1)
        if isinstance(edate, type(datetime.datetime)):
            edate = datetime.datetime(2018, 6, 1)
        dates = pd.date_range(sdate, edate)
    assert folder != type(None), 'ABORTING: folder for volcano files needed'
    # Write out the new files to the same folder if a new one is not given
    if isinstance(output_folder, type(None)):
        output_folder = folder
    # Helper functions for processing volcano files

    def get_volc_subfolder4dt(folder=None, dt=None):
        """ Get the volcano folder for a specific datetime """
        year = dt.year
        month = dt.month
        folder = '{}/{}/{:0>2}/'.format(folder, year, month)
        return folder

    def get_volc_filename4dt(dt):
        """ Get the filename for a specific datetime """
        dt_str = dt.strftime('%Y%m%d')
        return 'so2_volcanic_emissions_Carns.{}.rc'.format(dt_str)

    def get_volc_file_lines4date(dt, folder=None):
        """ Open a volcano file for given date and return its lines """
        # Get string contain filename string and another for folder
        filename = get_volc_filename4dt(dt=dt)
        folder2use = get_volc_subfolder4dt(dt=dt, folder=folder)
        # Extract lines
        lines = read_lines_from_txt_file(filename=filename,
                                         folder=folder2use)
        return lines

    def rm_eruption_from_volc_file_lines(lines, skiplines=4, dt_str='',
                                         verbose=True, debug=False):
        """
        Remove the eruptive volcano lines from the HEMCO file

        Notes
        ------
         - The columns (and units) in the files are:
            LAT (-90,90), LON (-180,180), SULFUR [kg S/s], ELEVATION [m],
            CLOUD_COLUMN_HEIGHT [m]
        """
        # Local variables
        pstr1 = "NOTE: rm'd emission of {} kg S/s ({}) - ({}N, {}E, {}m, {}m)"
        pstr2 = 'WARNING: line inc. without check: {}'
        NewLines = []
        # Now loop lines and only save those without eruptive volcanoes
        for n_line, line in enumerate(lines):
            include_line = True
            if (n_line+1) > skiplines:
                try:
                    # Get elevation and emission height
                    tmp_line = line.strip().split()
                    LAT = tmp_line[0]
                    LON = tmp_line[1]
                    S = tmp_line[2]  # SULFUR
                    ELEV = tmp_line[3]  # ELEVATION
                    CLOUD = tmp_line[4]  # CLOUD_COLUMN_HEIGHT
                    # If not equal, then EXCLUDE the line
                    if ELEV != CLOUD:
                        if verbose:
                            print(pstr1.format(S, dt_str, LAT, LON, ELEV,
                                               CLOUD))
                        include_line = False
                except IndexError:
                    if debug:
                        print(pstr2.format(line.strip()))
            else:
                if debug:
                    print(pstr2.format(line.strip()))
            if include_line:
                NewLines += [line]
        return NewLines

    # - Now loop dates and re write files
    for dt in dates:
        dt_str = dt.strftime('%Y/%m/%d')
        # Get lines of volcano file for date
        lines = get_volc_file_lines4date(dt=dt, folder=folder)
        # Remove the eruptive volcanoes from the file
        NewLines = rm_eruption_from_volc_file_lines(lines=lines, dt_str=dt_str)
        if len(lines) != len(NewLines):
            print('WARNING: # of lines updated for {}'.format(dt_str))
#        SubFolder='VOLCANO_NO_ERUPT/v2019-08/'
        # Get string contain filename string
        filename = get_volc_filename4dt(dt=dt)
        folder2use = get_volc_subfolder4dt(dt=dt, folder=output_folder)
        # Save out the new file
        write_lines2txt_file(NewLines, folder=folder2use, filename=filename)
