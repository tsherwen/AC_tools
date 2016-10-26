# Downlaod the example files from an external source

import os
from urllib2 import urlopen

example_dir = "data"
example_web = "http://atmosviz1.york.ac.uk/~bn506/data/AC_tools/"


file_list = [
"test.nc",
"tracerinfo.dat",
"diaginfo.dat",
]


if not os.path.exists(example_dir):
    os.makedirs(example_dir)

for _file in file_list:
    if _file == "test.nc":
        new_filename = os.path.join(example_dir, "ctm.nc")
    else:
        new_filename = os.path.join(example_dir, _file)

    if not os.path.isfile( new_filename ):
        new_file = open(new_filename, 'wb')
        file_url = example_web + _file
        print file_url
        file_data = urlopen( file_url ).read()
        new_file.write( file_data )
        new_file.close()
    
    

    


