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
# retain back compatibility for PyGChem
try:
  from pygchem import datasets
except:
  import pygchem.datafields as datasets

def convert_to_netCDF(folder='none',filename='ctm.nc',\
                         bpch_file_list=None, remake=False):
   """    
    Converts GEOS-Chem ctm.bpch output file(s) to NetCDF
   """   
   # Check if file already exists and warn about remaking
   from bpch2netCDF import get_folder
   folder = get_folder(folder)
   output_file = folder + '/' + filename

   # If the netCDf file already exists dont overwrite it without remake=True.
   if not remake:
       if os.path.exists(output_file):
           logging.warning(output_file + ' already exists. Not recreating.')
           return
       
   # By default look inside the folder for any files
   if bpch_file_list==None:
       bpch_files = glob.glob( folder + '/*.bpch*' )
       if len(bpch_files) == 0:
          bpch_files = glob.glob( folder + '/*trac_avg*' )
          if len(bpch_files) == 0:
               logging.error("No bpch files found in "+folder)
               raise IOError(folder + " contains no bpch files.")
   # Confirm the specified bpch files are there.
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
   print bpch_files
   bpch_data = datasets.load(bpch_files)

   # Save the netCDF file
#   iris.fileformats.netcdf.save(data, output_file)
   datasets.save( bpch_data, output_file )

   return

def get_folder(folder):
   """
    Get name of folder that contains ctm.bpch data from command line 
   """
   if folder=='none':
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

