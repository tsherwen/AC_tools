#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Functions for use with atmospheric observations (e.g. from FAAM BAE146 aircraft)  

Use help(<name of function>) to get details on a particular function.
"""
from bs4 import BeautifulSoup
import requests
import re
import wget
import urllib
import json



def get_FAAM_locations_as_df(flight_ID='C225'):
    """
    Retive the FAAM BAE146 position (current of historic) from the html website by flight ID
    """
    # What is the root URL for the data?
    URL = 'https://www.faam.ac.uk/gluxe/position/query?flight={}'.format(flight_ID)
    # Parse the URL via requests
    f = urllib.request.urlopen( URL )
    soup = BeautifulSoup( f )
    s = soup.get_text()
    # Parse the data a JSON string
    json_acceptable_string = s.replace("'", "\"")
    d = json.loads(json_acceptable_string)
    # Return as a dataframe
    return pd.DataFrame( d )
