#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Functions for use with atmospheric observations (e.g. from FAAM BAE146 aircraft)

Use help(<name of function>) to get details on a particular function.
"""
#from bs4 import BeautifulSoup
#import requests
import re
#import wget
import urllib
import json
try:
    import requests
except ImportError:
    print("WARNING: failed to import Python module 'requests'")

def get_FAAM_locations_as_df(flight_ID='C225'):
    """
    Retive the FAAM BAE146 position (current of historic) from the html website by flight ID
    """
    # What is the root URL for the data?
    URL = 'https://www.faam.ac.uk/gluxe/position/query?flight={}'.format(
        flight_ID)
    # Parse the URL via requests
    f = urllib.request.urlopen(URL)
    soup = BeautifulSoup(f)
    s = soup.get_text()
    # Parse the data a JSON string
    json_acceptable_string = s.replace("'", "\"")
    d = json.loads(json_acceptable_string)
    # Return as a dataframe
    return pd.DataFrame(d)


def sort_sites_by_lat(sites):
    """
    Order given list of GAW sties by latitudes
    """
    # Get info
    vars = [gaw_2_loc(s) for s in sites]  # lat, lon, alt, TZ
    # Sort by lat, index orginal sites list and return
    lats = [i[0] for i in vars]
    slats = sorted(lats)[::-1]
    return [sites[i] for i in [lats.index(ii) for ii in slats]]