#!/usr/bin/python
""" 
This script analysis a folder containing bpch files and outputs the results
in a single netCDF file in the folder.

This allows for significantly faster and easier input of data, 
and more common anaylsis techniques like pandas without extra 
post processing. 
"""

import logging
import sys
import glob
import os
import netCDF4
import iris
# retain back compatibility for PyGChem
try:
  from pygchem import datasets
except:
  import pygchem.datafields as datasets

def convert_to_netCDF(folder=None,filename='ctm.nc', bpch_file_list=None, \
        remake=False, hemco_file_list=None, verbose=True, \
        bpch_file_type="*.ctm.nc"):
    """
    Converts GEOS-Chem outputs to netCDF

    Parameters
    ----------
    folder (str): specify the folder you want to use - defaults to cwd
    filename (str):  specific the netCDF filename you want to use
    bpch_file_list (list): list the bpch files you want to use
    remake (boolean): Overwrite any old files (default=False)

    Notes
    -----
    Setup for:
    - bpch_to_netCDF
    - hemco_to_netCDF
    - planeflight_to_netCDF
    """
    logging.debug( "Convert to netCDF called with folder={},".format(folder)+ \
        " bpch_file_type={}/filename={}".format(bpch_file_type, filename) )

#    try:
    bpch_to_netCDF( folder=folder, filename=filename, 
            bpch_file_list=bpch_file_list, remake=remake, 
            file_type = bpch_file_type, verbose=verbose)
#    except:
#        logging.error("Could not convert bpch to netCDF in {_dir}"\
#                .format(_dir=folder))
#    try:
#        hemco_to_netCDF( folder, hemco_file_list, remake)
#    except:
#        logging.warning("Could not convert hemco to netCDF in {_dir}"\
#                .format(_dir=folder))

    return

def hemco_to_netCDF( folder, hemco_file_list=None, remake=False ):
    """
    Conbine HEMCO diagnostic output files to a single NetCDF file.

    Parameters
    ----------
    remake (boolean): overwrite existing NetCDF file

    """

    from bpch2netCDF import get_folder
    folder = get_folder(folder)
    output_file = os.path.join(folder, 'hemco.nc')

    # If the hemco netCDF file already exists then quit unless remake=True
    if not remake:
        if os.path.exists(output_file):
            logging.warning( output_file + ' already exists, not remaking')
            return

    logging.info("Combining hemco diagnostic files")

    # By default look for any files that look like hemco diagnostic files:
    # Look for all hemco netcdf files then remove the restart files.
    if hemco_file_list==None:
        hemco_files = glob.glob(folder + '/*HEMCO*.nc')
        for filename in hemco_files:
            if "restart" in filename:
                hemco_files.remove(filename)

    else:
      file_list = []
      for hemco_file in hemco_file_list:
         full_path = os.path.join(folder, hemco_file)
         if not os.path.exists(full_path):
            logging.error(full_path + " could not be found")
            raise IOError("{path} could not be found".format(path=full_path))
         file_list.append(full_path)
      hemco_files = file_list

    if len(hemco_files)==0:
        logging.warning("No hemco diagnostic files found in {_dir}"\
                .format(_dir=folder))
    else:
        logging.debug( "The following hemco files were found:")
        logging.debug( str(hemco_files) )

    # Use iris cubes to combine the data into an output file

    hemco_data = iris.load(hemco_files)
    # Concatanate the times.
    hemco_data = hemco_data.concatenate()
    iris.save( hemco_data, output_file)

    logging.info( str(hemco_data) )
    logging.info("Hecmo file created at {file}".format(file=output_file))
    return





    
 
def bpch_to_netCDF(folder=None, filename='ctm.nc', bpch_file_list=None, \
        remake=False, filetype="*ctm.bpch*", verbose=False, **kwargs):

   """    
   Converts GEOS-Chem ctm.bpch output file(s) to NetCDF

   Parameters
   ----------
   folder (str): working directory for data files
   filename (str): name to give created NetCDF
   bpch_file_list (list): list of files to convert 
   remake (boolean): overwrite existing NetCDF file
   filetype (str): string with wildcards to match filenames 
   ( e.g. *ctm.bpch*,*ts*bpch* )
   verbose (boolean): print (minor) logging to screen
   
   Returns
   -------
   (None) saves a NetCDF file to disk

   """   

   # Check if file already exists and warn about remaking
   from bpch2netCDF import get_folder
   folder = get_folder(folder)
   output_file = os.path.join(folder, filename)

   # If the netCDf file already exists dont overwrite it without remake=True.
   if not remake:
       if os.path.exists(output_file):
           logging.warning(output_file + ' already exists. Not recreating.')
           return
       
   # Look for files if file list is not provided.
   if isinstance( bpch_file_list, type(None) ):
       logging.debug("Searching for the following bpch filetype: {filetype}"\
                .format(filetype=filetype))
       bpch_files = glob.glob( folder + '/' + filetype )
       if len(bpch_files) == 0:
           logging.error("No bpch files found in "+folder)
           raise IOError(folder + " contains no bpch files.")

   # use the specified files.
   else:
       file_list = []
       for bpch_file in bpch_file_list:
         full_path = folder + '/' + bpch_file
         if not os.path.exists(full_path):
            logging.error(full_path + " could not be found")
            raise IOError("Full path could not be found")
         file_list.append(full_path)
       bpch_files = file_list

   # Open the bpch files
   logging.debug( "The following bpch files were found:")
   logging.debug( str(bpch_files) )
   if verbose:
        print "Creating a netCDF from {} file(s).".format(len(bpch_files))+\
            "This can take some time..."
   bpch_data = datasets.load(bpch_files)

   # Save the netCDF file
#   iris.fileformats.netcdf.save(data, output_file)
   datasets.save( bpch_data, output_file )
   logging.info( "A netCDF file has been created with the name {ctm}".format(ctm=output_file)) 
   return

def get_folder(folder):
   """
    Get name of folder that contains ctm.bpch data from command line 
   """
   if isinstance( folder, type(None) ):
      # getting the folder location from system argument
      if len(sys.argv)<=1:
         logging.warning( "No folder location specified for the data")
         folder = os.getcwd()
      else:
         folder = str(sys.argv[1])

   # Check folder exists
   if not os.path.exists( folder ):
      print "Folder does not exist"
      print folder
      sys.exit()


   return folder;
   
if __name__ == "__main__":
   convert_to_netCDF()
   print "Complete"

