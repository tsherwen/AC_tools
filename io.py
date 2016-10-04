#!/usr/bin/python
"""
This script analysis a folder containing bpch files and outputs the results
in a single netCDF file in the folder.

This allows for significantly faster and easier input of data, 
and more common anaylsis techniques like pandas without extra 
post processing. 
"""

def open_netCDF(folder='none',filename='ctm.nc', bpch_file_names=None,
                remake=False):
    """
    Opens the netCDF file from a GEOS-Chem run.
    Converts all .bpch files in a folder to a netCDF file if required.
    The Default folder is the current folder.
    The default ouptut filename is ctm.nc.
    Returns a netCDF4 Dataset object.
    """

    import logging
    import os
    if (folder == 'none'):
        folder = os.getcwd()
    # Strip trailing "/"
    elif (folder[-1]=="/"):
        folder = folder[:-1]
        
    netCDF_filename = os.path.join(folder, filename)
    
    if remake:
        try:
            os.remove(netCDF_filename)
            logging.info("netCDF file deleted for remake")
        except:
            logging.info("Tried remaking but could not remove old file")


           

    logging.info('Opening netCDF file: ' + netCDF_filename)


    from netCDF4 import Dataset
    # Try to open a netCDF file
    try:
        # Hope for success
        netCDF_data = Dataset(netCDF_filename)
        logging.info("netCDF file opened successfuly")
    except:
        # If no netCDF file loaded, try making one from bpch
        logging.debug('No netCDF file found. Attempting to create one.')
    
        # Confirm aux files exist (tracerinfo.dat)
        # If they are not then the netCDF variable names can be bugged.
        # Code to do this check goes here...


        try:
            from pygchem import datasets
        except:
            import pygchem.datafields as datasets
        import glob
        if (bpch_file_names==None):
            bpch_files = glob.glob(folder + '/*.bpch' )
        else:
            bpch_files = []
            for bpch_file in bpch_file_names:
                bpch_location = folder + '/' + bpch_file
                if not os.path.isfile(bpch_location):
                    raise IOError('{_file} not found'.format(_file=bpch_location))
                bpch_files.append(bpch_location)
        logging.debug('Found ' + str(len(bpch_files)) + ' bpch files.')
        if (len(bpch_files) == 0):
            logging.error('No bpch files found.')
            raise IOError('Cannot find bpch files in {folder}'\
                        .format(folder=folder))
            return
        pygchem_data = datasets.load(bpch_files)
        # Convert the dataset to an iris cube and export as netCDF
        logging.debug('Creating a netCDF file from bpch.')
        print "Creating a netCDF file from bpch"
        datasets.save( pygchem_data, netCDF_filename)
        # Save the dataset to disk


        # Open the netCDF dataset.
        netCDF_data = Dataset( netCDF_filename )
        
        if (netCDF_data == None):
            logging.error('Error creating netCDF file from bpch.')
            raise IOError('Error creating netCDF file from bpch.')
    
        # Confirm that the netCDF file time size is the same size 
        #  as the size of the list
        # To-do: There is probably a cleaner way to do this.
        if not (len(bpch_files) == 1):
            if not (len(netCDF_data.variables['time']) == len(bpch_files)):
                logging.error('Incorrect amount of timesets from bpch files.'\
                              'Could be due to incomplete bpch files')
                # Find any potential small files
                for bpch_file in bpch_files:
                    if (os.path.getsize( bpch_file ) < 1000):
                        logging.info('{filename} looks rather small...'\
                                    .format(filename=bpch_file) )
        #        raise IOError('Incorrect amount of timesets from bpch files.')
                            
                                
            


    return netCDF_data
    
