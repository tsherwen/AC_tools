# Downlaod the example files from an external source

import os
from urllib2 import urlopen
import logging

data_dir = "../data"
data_url = "http://atmosviz1.york.ac.uk/~bn506/data/AC_tools/"


file_list = [
"HEMCO_Diagnostics.nc",
"ctm.nc",
"diaginfo.dat",
"test.bpch",
"test.log",
"tracerinfo.dat",
"LANDMAP/LANDMAP_LWI_ctm/ctm.nc",
"LANDMAP/LANDMAP_LWI_ctm_025x03125/ctm.nc",
"LANDMAP/LANDMAP_LWI_ctm_05x0666/ctm.nc",
"LANDMAP/LANDMAP_LWI_ctm_2x25.lod/ctm.nc",
"LANDMAP/LANDMAP_LWI_ctm_2x25/ctm.nc",
"LANDMAP/LANDMAP_ctm_2x25/ctm.nc",
]


if not os.path.exists(data_dir):
    os.makedirs(data_dir)

for _file in file_list:
#    if _file == "test.nc":
#        new_filename = os.path.join(data_dir, "ctm.nc")
#    else:
    new_filename = os.path.join(data_dir, _file)

    if not os.path.isfile( new_filename ):
        try:
            os.makedirs(os.path.dirname(new_filename))
        except:
            logging.error("Could not create folder for {file}".format(file=new_filename))
        logging.debug( new_filename + " not found. Downloading now.")
        new_file = open(new_filename, 'wb')
        file_url = data_url + _file
        try:
            logging.debug("downloading from {url}".format(url=file_url))
            file_data = urlopen( file_url ).read()
            new_file.write( file_data )
            logging.debug( new_filename+" downloaded.")
        except:
            logging.error("Failed to download {url}".format(url=file_url))
        new_file.close()
    
    

    


