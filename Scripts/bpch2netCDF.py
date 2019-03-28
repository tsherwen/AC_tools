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
try:
    import iris
except ImportError:
    print('WARNING iris not imported')
# retain back compatibility for PyGChem
try:
    if (sys.version_info.major <= 2):
        import pygchem
        if pygchem.__version__ == '0.2.0':
            import pygchem.diagnostics as gdiag
        else:
            try:
                from pygchem import datasets
            except:
                import pygchem.datafields as datasets
except ImportError:
    print('pygchem not imported!')


def convert_to_netCDF(folder=None, filename='ctm.nc', bpch_file_list=None,
                      remake=False, hemco_file_list=None, verbose=True,
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
    logging.debug("Convert to netCDF called with folder={},".format(folder) +
                  " bpch_file_type={}/filename={}".format(bpch_file_type, filename))

#    try:
    bpch_to_netCDF(folder=folder, filename=filename,
                   bpch_file_list=bpch_file_list, remake=remake,
                   file_type=bpch_file_type, verbose=verbose)
#    except:
#        logging.error("Could not convert bpch to netCDF in {_dir}"\
#                .format(_dir=folder))
#    try:
#        hemco_to_netCDF( folder, hemco_file_list, remake)
#    except:
#        logging.warning("Could not convert hemco to netCDF in {_dir}"\
#                .format(_dir=folder))

    return


def hemco_to_netCDF(folder, hemco_file_list=None, remake=False):
    """
    Conbine HEMCO diagnostic output files to a single NetCDF file.

    Parameters
    ----------
    remake (boolean): overwrite existing NetCDF file

    """
    if __package__ is None:
        from .bpch2netCDF import get_folder
    else:
        from .bpch2netCDF import get_folder
    folder = get_folder(folder)
    output_file = os.path.join(folder, 'hemco.nc')

    # If the hemco netCDF file already exists then quit unless remake=True
    if not remake:
        if os.path.exists(output_file):
            logging.warning(output_file + ' already exists, not remaking')
            return

    logging.info("Combining hemco diagnostic files")

    # By default look for any files that look like hemco diagnostic files:
    # Look for all hemco netcdf files then remove the restart files.
    if hemco_file_list == None:
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
                raise IOError(
                    "{path} could not be found".format(path=full_path))
            file_list.append(full_path)
        hemco_files = file_list

    if len(hemco_files) == 0:
        logging.warning("No hemco diagnostic files found in {_dir}"
                        .format(_dir=folder))
    else:
        logging.debug("The following hemco files were found:")
        logging.debug(str(hemco_files))

    # Use iris cubes to combine the data into an output file

    hemco_data = iris.load(hemco_files)
    # Concatanate the times.
    hemco_data = hemco_data.concatenate()
    iris.save(hemco_data, output_file)

    logging.info(str(hemco_data))
    logging.info("Hecmo file created at {file}".format(file=output_file))
    return


def bpch_to_netCDF(folder=None, filename='ctm.nc', bpch_file_list=None,
                   remake=False, filetype="*ctm.bpch*",
                   check4_trac_avg_if_no_ctm_bpch=True, backend='PyGChem',
                   verbose=False, **kwargs):
    """
    Converts GEOS-Chem ctm.bpch output file(s) to NetCDF

    Parameters
    ----------
    folder (str): working directory for data files
    filename (str): name to give created NetCDF
    bpch_file_list (list): list of files to convert
    remake (boolean): overwrite existing NetCDF file
    filetype (str): string with wildcards to match filenames
    ( e.g. *ctm.bpch*, trac_avg.*, or *ts*bpch* )
    verbose (boolean): print (minor) logging to screen

    Returns
    -------
    (None) saves a NetCDF file to disk
    """
    import os
    # Check if file already exists and warn about remaking
    if __package__ is None:
        from .bpch2netCDF import get_folder
    else:
        from .bpch2netCDF import get_folder
    folder = get_folder(folder)
    output_file = os.path.join(folder, filename)

    # If the netCDf file already exists dont overwrite it without remake=True.
    if not remake:
        if os.path.exists(output_file):
            logging.warning(output_file + ' already exists. Not recreating.')
            return

    # Look for files if file list is not provided.
    if isinstance(bpch_file_list, type(None)):
        logging.debug("Searching for the following bpch filetype: {filetype}"
                      .format(filetype=filetype))
        bpch_files = glob.glob(folder + '/' + filetype)
        # Also check if directory contains *trac_avg* files, if no ctm.bpch
        if (len(bpch_files) == 0) and check4_trac_avg_if_no_ctm_bpch:
            filetype = '*trac_avg*'
            logging.info('WARNING! - now trying filetype={}'.format(filetype))
            bpch_files = glob.glob(folder + '/' + filetype)
        # Raise error if no files matching filetype
        if len(bpch_files) == 0:
            logging.error("No bpch files ({}) found in {}".format(filetype,
                                                                  folder))
            raise IOError("{} contains no bpch files.".format(folder))

    # Use the specified files.
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
    logging.debug("The following bpch files were found (n={}):"
                  .format(len(bpch_files)))
    logging.debug(str(bpch_files))
    if verbose:
        print(("Creating a netCDF from {} file(s).".format(len(bpch_files)) +
               " This can take some time..."))
    if backend == 'PyGChem':
        # Load all the files into memory
        bpch_data = datasets.load(bpch_files)
        # Save the netCDF file
        datasets.save(bpch_data, output_file)
    elif backend == 'xbpch':
        import xbpch
        # Load all the files into memory (as xarray dataset object)
        ds = xbpch.open_mfbpchdataset(bpch_files)
        # save through xarray dataset object
        ds.to_netcdf(folder+filename,
                        unlimited_dims={'time_counter': True})
    elif backend == 'iris':
        #    iris.fileformats.netcdf.save(data, output_file)
        print('WARNING NetCDF made by iris is non CF-compliant')
    elif backend == 'PNC':
        import PseudoNetCDF as pnc
        import xarray as xr
        if len(bpch_files) == 1:
            bpch_to_netCDF_via_PNC(filename=filename,
                                   output_file=output_file, bpch_file=bpch_files[0])
        # Individually convert bpch files if more than one file
        if len(bpch_files) > 1:
            for n_bpch_file, bpch_file in enumerate(bpch_files):
                bpch_to_netCDF_via_PNC(filename=filename,
                                       output_file='TEMP_{}_'.format(
                                           n_bpch_file)+filename,
                                       bpch_file=bpch_file)
            # - Combine the NetCDF files with xarray
            TEMP_ncfiles = glob.glob(folder+'TEMP_*_'+filename)
            # Open files with xarray
            ds_l = [xr.open_dataset(i) for i in TEMP_ncfiles]
            # Make sure the time dimension is unlimitetd
            ds = xr.concat(ds_l, dim='time')
            # Now save the combined file
            ds.to_netcdf(folder+filename,
                         unlimited_dims={'time_counter': True})
            # Remove the temporary files
            for TEMP_ncfile in TEMP_ncfiles:
                os.remove(TEMP_ncfile)

    logging.info("A netCDF file has been created with the name {ctm}"
                 .format(ctm=output_file))
    return


def bpch_to_netCDF_via_PNC(format='bpch2', filename='ctm.nc',
                           output_file=None, bpch_file=None, folder=None):
    """ Convert bpch to NetCDF using PNC as backend """
    import PseudoNetCDF as pnc
    # Load the file into memory
    infile = pnc.pncopen(bpch_file, format=format)
    # Kludge - reduce DXYP_DXYP dims online
    dxyp = infile.variables['DXYP_DXYP']
    # Surface area should have time dim, if fit does remove it.
    if len(dxyp.shape) == 4:
        dxyp.dimensions = dxyp.dimensions[1:]
        infile.variables['DXYP_DXYP'] = dxyp
    # Now write file to disc
#    pnc.pncwrite(infile, folder+filename)
    pnc.pncwrite(infile, output_file)


def get_folder(folder):
    """
    Get name of folder that contains ctm.bpch data from command line
    """
    if isinstance(folder, type(None)):
       # getting the folder location from system argument
        if len(sys.argv) <= 1:
            logging.warning("No folder location specified for the data")
            folder = os.getcwd()
        else:
            folder = str(sys.argv[1])

    # Check folder exists
    if not os.path.exists(folder):
        print("Folder does not exist")
        print(folder)
        sys.exit()

    return folder


if __name__ == "__main__":
    convert_to_netCDF()
    print("Complete")
