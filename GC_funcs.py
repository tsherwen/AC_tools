"""
Functions for GeosChem netCDF files.
Based on TS AC_tools.funcs4GEOSC
All opperations are performed on a netCDF file created by AC_tools.io
"""

import logging
import numpy as np

def get_variable_data( netCDF_file, variable ):
    """
    Get variable data from netCDF file.
    """
    logging.info('Getting {variable} data from netCDF file'.format(variable=variable))

    try:
        variable_data = netCDF_file.variables[variable]
        dimensions = variable_data.dimensions
        logging.info('Found data with dimensions {dimensions}'.format(dimensions=dimensions))
    except:
        logging.error('Could not find variable {variable}'.format(variable=variable))
        raise IOError('Could not find variable {variable}'.format(variable=variable))
    return variable_data

def get_surface_area( netCDF_file ):
    """
    Optains a 2D map of gridBOx surfice areas.
    """
    from .GC_funcs import get_variable_data
    logging.info( 'Getting surfice area from netCDF file' )

    surface_area = get_variable_data( netCDF_file, 'DXYP__DXYP' )

    if not (len(surface_area.shape) == 2):
        logging.error('Problem obtaining the surfice area')

    return surface_area;

def get_variables( netCDF_file, show=True ):
    """
    Print a list of variables included in the netCDF file.
    """
    logging.info("Getting the list of netCDF variables")

    variables = []
    for variable in netCDF_file.variables:
        variables.append(variable)
    if show==True:
        for variable in variables:
            print variable
    return variables;


def get_land_map( netCDF_file ):
    """
    Get the land, water, ice (LWI) from GEOS-Chem with integers for land(1) and water(0).
    Ice fraction is provided as a fractional vaule.
    """
#    from .GC_funcs import get_variable_data
#    logging.info('Getting the land map.')
#    land_map = get_variable_data(netCDF_file, 'LANDMAP__LWI')
#    return land_map
    return

def get_air_mass( netCDF_file ):
    """
    Get the air mass of every gridbox.
    """

    from .GC_funcs import get_variable_data
    logging.info('Getting the air mass.')
    air_mass = get_variable_data(netCDF_file, 'BXHGHT_S__AD')

    return air_mass
    

def get_trop_time( netCDF_file ):
    """
    Get the percentage of time each gridbox spends in the troposphere
    """

    from .GC_funcs import get_variable_data
    trop_time_variable_name = "TIME_TPS__TIMETROP"
    trop_time = get_variable_data( netCDF_file, trop_time_variable_name )
    return trop_time

def get_species_data( netCDF_file, species ):
    """
    Gets the species data in mixing ratio from the netCDF file.
    """
    logging.info( 'Getting {species} data.'.format(species=species) )

    assert isinstance(species, str), 'Species provided is not a string'
    from .GC_funcs import get_variable_data
    variable_name = 'IJ_AVG_S__' + species
    variable_data = get_variable_data( netCDF_file, variable_name)
    return variable_data

def get_species_rmm( species_name ):
    """
    Gets a species rmm.
    """
    from .species import species as get_species
    
    species = get_species( species_name)
    try:
        species_rmm = species.RMM
    except:
        logging.warning("No RMM found for {species_name}. Assuming mass of 1.0"\
            .format(species_name=species_name))
        species_rmm = 1.0
    return species_rmm
    
    

def get_tropospheric_species_mass( netCDF_file, species ):
    """
    Return the mass of a species in each gridbox in gramms.
    """
    
    from .GC_funcs import get_species_data
    from .GC_funcs import get_air_mols
    from .GC_funcs import get_trop_time
    from .GC_funcs import get_species_rmm

    logging.info('Getting the tropospheric {species} mass.'.format(species=species))

    species_data = get_species_data( netCDF_file, species )
    trop_time = get_trop_time( netCDF_file )
    species_rmm = get_species_rmm(species)

    air_moles = get_air_mols( netCDF_file)
    species_moles = np.multiply(air_moles, np.divide(species_data,1E9))
    species_mass = np.multiply( species_moles, species_rmm )

    if len(species_mass.shape)==3:
        trop_species_mass = np.multiply( species_mass[:,:,:38], trop_time )
    elif len(species_mass.shape)==4:
        trop_species_mass = np.multiply( species_mass[:,:,:,:38], trop_time )
    else:
        logging.error("species_mass is the wrong shape")
    return trop_species_mass


def get_air_mols( netCDF_file ):
    """"
    Get the number of molecules in each gridbox.
    """
    from .GC_funcs import get_air_mass

    air_mass = get_air_mass( netCDF_file )
    # Air density for N2 + O2
    air_moles = np.divide(np.multiply(air_mass,1E3) , (0.78*28.0 + 0.22*32.0))

    return air_moles


def get_tropospheric_burden( netCDF_file, variable ):
    """
    Get the tropospheic burden in Tg from the provided netCDf file.
    """
    from .GC_funcs import get_tropospheric_species_mass

    trop_species_mass = get_tropospheric_species_mass( netCDF_file, variable )

    # Average over all times if time dimension exists
    if len(trop_species_mass.shape)==4:
        trop_species_mass = np.mean(trop_species_mass, axis=0)


    # Get the total in Tg.
    total_trop_species_mass = np.sum(trop_species_mass) / 1E12

    return float(total_trop_species_mass)

def get_volume( netCDF_file ):
    """
    Get the volume of ell gridboxes in m^3.
    """
    from .GC_funcs import get_variable_data

    height = get_variable_data(netCDF_file, "BXHGHT_S__BXHEIGHT")
    area = get_variable_data(netCDF_file, "DXYP__DXYP")

    # Move time to the first dimension for easier multiplication.
    height = np.rollaxis(height[:],-1,0)
    volume = np.multiply( height, area )
    # Move time back to final axis
    volume = np.rollaxis(volume, 0, np.size(volume.shape))    

    if (len(volume.shape) == 5):
        volume = volume[:,:,:,:,0]

    return volume

def get_total_OH_PL(netCDF_file):
    """
    get the total amount of tropospheric OH production and loss.
    """
    from .GC_funcs import get_variable_data
    # [molec / cm3]
    OH_PL = get_variable_data( netCDF_file, 'CHEM_L_S__OH')

    from .GC_funcs import molec_cm3_to_Tg
    from .GC_funcs import get_volume

    # [m3]
    volume = get_volume(netCDF_file)

    OH_RMM = 17.0
    # [Tg]
    arr = molec_cm3_to_Tg( OH_PL, volume, OH_RMM )

    from .GC_funcs import get_trop_time
    trop_time = get_trop_time( netCDF_file )

    # [Tg]
    arr = np.multiply(arr, trop_time)

    if len(arr.shape)==4:
        arr = np.mean(arr, axis=0)
    arr = np.sum(arr)

    return arr


    
    



#def get_total_HO2_PL(netCDF_file):

    
def get_tropospheric_PL(netCDF_file, group_name, group_RMM):
    """
    Get the amount per gridbox of tropsopheric production for 
    a group in Tg per gridbox.
    """
    logging.info("Getting the tropospheric prod/loss for {var}"\
                .format(var=group_name))

    from .GC_funcs import get_variable_data
    from .GC_funcs import get_air_mass
    from .GC_funcs import get_trop_time
    from .GC_funcs import get_air_mols
    from .GC_funcs import get_volume
    from .GC_funcs import molec_cm3_s_to_Tg_year

    assert isinstance(group_name, str), "Invalid PL group name. Not a string."
    assert isinstance(group_RMM, float), "Group RMM not a float."

    variable_name = "PORL_L_S__"+group_name
    
    arr = get_variable_data( netCDF_file, variable_name )
    air_mols = get_air_mols( netCDF_file )
    air_mass = get_air_mass( netCDF_file )
    trop_time = get_trop_time( netCDF_file )
    volume = get_volume( netCDF_file ) # Get in cm^3


    arr = molec_cm3_s_to_Tg_year( arr, volume, group_RMM )



    # Get only the tropospheric part
    arr = np.multiply( arr, trop_time )
    
    return  arr

def molec_cm3_to_Tg( arr, volume, group_RMM):
    """
    converts an array from molec / cm3 to Tg .
    arr = array_input (molec / cm^3 )
    volume = volume of gridboxs (m^3)
    group_RMM = relative molecular mass of the species ( # )
    """

    assert (len(arr.shape) in [3,4]),\
         "The array shape of {shape} is not supported"\
         .format(shape=arr.shape)

    assert (len(volume.shape) in [3,4]),\
         "The volume shape of {shape} is not supported"\
         .format(shape=volume.shape)
    assert isinstance( group_RMM, float), "Group RMM is not a float"
    
    # Convert volume [m3] to [cm3]
    volume = np.multiply( volume, 1E6 )
    
    # Convert from molec/cm3 to molec
    if len(volume.shape)==3:
        arr = np.multiply( arr, volume[:,:,:38])
    elif len(volume.shape)==4:
        arr = np.multiply( arr, volume[:,:,:,:38])


    # Convert from molec to moles
    avagadro = 6.0221409e+23
    arr = np.divide( arr, avagadro )

    # Convert from moles to gram
    arr = np.multiply(arr, group_RMM) 

    # Convert from gram to TG
    arr = np.divide( arr, 1E12 )

    return arr

def molec_cm3_s_to_Tg_year( arr, volume, group_RMM ):
    """
    converts an array from molec / cm3 / s to Tg year.
    arr = array
    volume = volume of gridboxs [m3]
    group_RMM = relative molecular mass of the species
    """

    from .GC_funcs import molec_cm3_to_Tg
    arr = molec_cm3_to_Tg( arr, volume, group_RMM )

    # Convert from gram/s to gram/y
    year = 365.25*24*60*60
    arr = np.multiply(arr, year)

    return arr
    

def get_tropospheric_total_PL(netCDF_file, group_name, group_RMM):
    """
    Get the total amount of tropospheric producution for a group in Tg.
    """
    logging.debug("Getting the total tropsopheric PL for {variable}".format(variable=group_name))

    from .GC_funcs import get_tropospheric_PL
    PL = get_tropospheric_PL( netCDF_file, group_name, group_RMM)
    if (len(PL.shape)==4):
        PL = np.mean(PL, axis=0)
    PL = np.sum(PL)

    PL = float(PL)


    return PL
    
def get_drydep(netCDF_file, species):
    """
    get the dry deposition rate in g/s
    """
    logging.debug("Getting the dry deposition for {variable}"\
                .format(variable=species))

    from .GC_funcs import get_variable_data
    from .GC_funcs import get_species_rmm
    
    variable_name = "DRYD_FLX__{variable}df".format(variable=species)

    # Get surface area [m2]
    area = get_variable_data(netCDF_file, "DXYP__DXYP")

    # Get molec/cm2/s
    drydep = get_variable_data( netCDF_file, variable_name )
    
    # get molec/m2/s
    drydep = np.multiply(drydep, 1E4)

    # get molec/s
    drydep = np.multiply( area, drydep )

    # get moles/s
    avagadro = 6.0221409e+23                                                    
    drydep = np.divide( drydep, avagadro )    

    # get g/s
    species_rmm = get_species_rmm(species)
    drydep = np.multiply( species_rmm, drydep)

    return drydep 

def get_annual_drydep( netCDF_file, group_name):
    """
    get the annual dry deposition rate in Tg/year,
    """
    logging.debug("Getting the annual dry deposition.")

    from .GC_funcs import get_drydep

    year = 365.25*24*60*60

    drydep = get_drydep( netCDF_file, group_name)
    # If time then average over time
    if len(drydep.shape)==3:
        drydep = np.mean(drydep, axis=0)#
    total_drydep = np.sum(drydep)

    # convert from per second to per year
    total_drydep = np.multiply(total_drydep, year)
    # convert from g to Tg
    total_drydep = np.divide(total_drydep , 1E12)

    return float(total_drydep)


    


