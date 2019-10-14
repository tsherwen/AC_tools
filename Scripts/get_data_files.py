"""
Download the example/reference data files for AC_tools from an external source
"""

import os
import logging
# temporarily restore urllib2 (as default urllib not working)
#import urllib
#import urllib2
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError
#except ModuleNotFoundError:
except ImportError:
    from urllib.request import urlopen
    from urllib.error import HTTPError

# Set up logging only if running as a script
if __name__ == '__main__':
    FORMAT = "%(levelname)8s - %(message)s   @---> %(filename)s:%(lineno)s  %(funcName)s()"
    logging.basicConfig(filename='AC_tools.log', filemode='w', level=logging.DEBUG,
                        format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)

# Location of files on York visible apache server (atmosviz1)
data_dir = "../data"
data_url = "https://webfiles.york.ac.uk/WACL/SHERWEN_TOMAS/AC_tools/"

# List of files to fetch
file_list = [
    "HEMCO_Diagnostics.nc",
    "ctm.nc",
    "diaginfo.dat",
    "test.bpch",
    "test.log",
    "tracerinfo.dat",
    "LM/LANDMAP_LWI_ctm_4x5/ctm.nc",
    "LM/LANDMAP_LWI_ctm_025x03125/ctm.nc",
    "LM/LANDMAP_LWI_ctm_025x03125_CH/ctm.nc",
    "LM/LANDMAP_LWI_ctm_025x03125_WA/ctm.nc",
    "LM/LANDMAP_LWI_ctm_05x0666/ctm.nc",
    "LM/LANDMAP_LWI_ctm_2x25/ctm.nc",
    "LM/LANDMAP_LWI_ctm_0125x0125/ctm.nc",
    "GEOS_ChemSpecies_fullchem_v0.1.0.csv"
]


def main():
    """
    Driver to download data/reference files for AC_tools
    """
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    for _file in file_list:
        new_filename = os.path.join(data_dir, _file)
        file_url = data_url + _file

        # If file does not exist donwload the file
        if not os.path.isfile(new_filename):
            logging.debug(new_filename + " not found. Downloading now.")
            download_file(new_filename, file_url)

        else:
            # If file exists make sure it is the correct size
            url_size = int(urlopen(file_url).info()['Content-Length'])
            file_size = int(os.stat(new_filename).st_size)
            if not url_size == file_size:
                logging.warning("{fn} appears to be the wrong size\
                        {size1} vs {size2}".format(
                    fn=new_filename, size1=url_size, size2=file_size))
                logging.warning("Redownloading now")
                download_file(new_filename, file_url)


def download_file(new_filename, file_url):
    """
    Download a file from a given URL
    """
    if not os.path.exists(os.path.dirname(new_filename)):
        try:
            os.makedirs(os.path.dirname(new_filename))
        except:
            logging.error(
                "Could not create folder for {file}".format(file=new_filename))

    try:
        new_file = open(new_filename, 'wb')
        logging.debug("downloading from {url}".format(url=file_url))
        prt_str = "Downloading file ({}), which may take some time."
        print( prt_str.format(new_filename) )
        file_data = urlopen(file_url).read()
        new_file.write(file_data)
        print("Download complete.")
        logging.debug(new_filename+" downloaded.")
        new_file.close()
    except HTTPError as error_code:
        if error_code.code == 404:
            logging.error("{The following was not found on the server:")
            logging.error("{url}".format(url=file_url))
        else:
            logging.error("Failed to get {url} with HTTP error {error_code}"
                          .format(url=file_url, error_code=error_code))
    except:
        logging.error("Failed to download {url}".format(url=file_url))


# Run whether as script or as import
main()
