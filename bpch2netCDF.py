#!/usr/bin/python
""" This script analysis a folder containing bpch files and outputs the results
        in a single netCDF file in the folder.

         This allows for significantly faster and easier input of data, 
         and more common anaylsis techniques like pandas without extra 
         post processing. """

def convert_to_netCDF(folder='none',filename='ctm.nc'):

   from bpch2netCDF import get_folder
   folder = get_folder(folder)

   # Get ctm.bpch ( inc. if named *trac_avg* )
   import glob
   bpch_files = glob.glob( folder + '/*ctm.bpch*' )
   if len(bpch_files) == 0:
      bpch_files = glob.glob( folder + '/*trac_avg*' )

   # Open the bpch files
   try:
      from pygchem import datasets
   except:
      import pygchem.datafields as datasets
   data = datasets.load(bpch_files)


   # Save the netCDF file
   output_file = folder + '/' + filename
   import iris
   iris.fileformats.netcdf.save(data, output_file)

   return

def get_folder(folder):

   import sys
   if folder=='none':
      # getting the folder location from system argument
      if len(sys.argv)<=1:
         print "No folder location specified for the data"
         sys.exit()
      else:
         folder = str(sys.argv[1])

   # Check folder exists
   import os
   if not os.path.exists( folder ):
      print "Folder does not exist"
      print folder
      sys.exit()


   return folder;
   
if __name__ == "__main__":
   convert_to_netCDF()
   print "Complete"

