#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic functions for use with GEOS-Chem/Data Analysis.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""
# - Required modules:
# I/O / Low level
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
# time
import time
import datetime as datetime
# math
from math import radians, sin, cos, asin, sqrt, pi, atan2
# geopands and rasterio are not back compatibly with python 2, so only try and import
try:
    import geopandas
except ImportError:
    print('WARNING: failed to import geopandas')
try:
    from rasterio import features
except ImportError:
    print('WARNING: failed to import rasterio')
from affine import Affine

# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from . core import *
from . variables import *


def get_country_mask(country='South Africa', res='2x2.5'):
    """
    Get a mask (1s) for a given country at the requested resolution
    """
    # Get AC_tools location, then set example data folder location
    import os
    import xarray as xr
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path = os.path.dirname(os.path.abspath(filename))
    folder = path+'/data/LM/LANDMAP_LWI_ctm_0125x0125/'
    # Get coords from LWI 0.125x0.125 data and remove the time dimension
    ds = xr.open_dataset(folder+'ctm.nc')
    ds = ds.mean(dim='time')
    # Add a raster mask for a country
    ds = add_raster_of_country2ds(ds, test_plot=True, country=country)
    # Only include states in the assignment
    ds = ds[['states']]
    # rrgrid to coarser resolution (e.g. 2x2.5)
    ds = regrid2coarse_res(ds, res=res)
    return ds


def add_raster_of_oceans2ds(ds, featurecla='ocean', set_all_regions2one=False,
                            test_plot=False, dpi=320):
    """
    Add raster outline of country to spatial dataset
    """
    # Get shapes for country
    shapes = get_shapes4oceans(featurecla=featurecla)
    # Add country's states as a layer
    ds[featurecla] = rasterize(shapes, ds.coords)
    # Test plot of this?
    if test_plot:
        from . plotting import quick_map_plot
        savename = 'spatial_plot_of_shapes4oceans_{}'.format(featurecla)
        quick_map_plot(ds, var2plot=featurecla, savename=savename)
    # set all the regions (e.g. counties/states) in a country to 1
    if set_all_regions2one:
        arr = ds[featurecla].values
        arr[np.where(~np.isnan(arr))] = 1
        ds[featurecla].values = arr
    return ds


def get_shapes4oceans(featurecla='ocean', rtn_group=False):
    """
    Get shapes (polygons) for oceans from Natural earth

    NOTES
    -------
     - data credit: NaturalEarth     https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_geography_marine_polys.zip
    """
    # location of data
    URL = "http://www.naturalearthdata.com/downloads/10m-physical-labels/"
    URL += "/10m-ocean/"
    # Shapefiles locally?
    # TODO - update to download automatically and store in AC_tools' data directory
#    shapefiles = 'ne_10m_ocean'
    shapefiles = 'ne_10m_geography_marine_polys'
    folder = '/mnt/lustre/users/ts551/labbook/Python_progs/'
    folder += '/AC_tools/data/shapefiles/{}'.format(shapefiles, shapefiles)
    group = geopandas.read_file(folder)
    # Just select state of interest
    choosen_group = group.query("featurecla == '{}'".format(featurecla))
    choosen_group = choosen_group.reset_index(drop=True)
    # Get the shapes
    shapes = zip(choosen_group.geometry, range(len(choosen_group)))
    if rtn_group:
        return choosen_group
    else:
        return shapes


def add_loc_ocean2df(df=None, LatVar='lat', LonVar='lon'):
    """
    Add the ocean of a location to dataframe
    """
    from geopandas.tools import sjoin
    # Get the shapes for the ocean
    featurecla = 'ocean'
    group = get_shapes4oceans(rtn_group=True, featurecla=featurecla)
    # Turn the dataframe into a geopandas dataframe
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df[LonVar], df[LatVar]))
    # Work out if any of the points are within the polygons
    pointInPolys = sjoin(gdf, group, how='left')
    # Check how many were assigned to a region
    Nnew = float(pointInPolys['name'].dropna().shape[0])
    N = float(df.shape[0])
    if N != Nnew:
        pstr = 'WARNING: Only {:.2f}% assigned ({} of {})'
        print(pstr.format((Nnew/N)*100, int(Nnew), int(N)))
    # Add the ocean assingnment back into the orginal dataframe
    df[featurecla] = pointInPolys['name'].values
    return df


def regrid2coarse_res(dsA, res='2x2.5', method='bilinear'):
    """
    Regrid a high resolution dataset to a lower resolution (e.g. 2x2.5)
    """
    import xesmf as xe
    # Hard code this for now
    grid2use = {
        '2x2.5': '2x2.5_deg_centre_GEOSChem',
        '4x5': '4x5_deg_centre_GEOSChem',
    }[res]
    # Get dictionary of grid coordinates
    grids = grids4reses()
    # Create a dataset to re-grid into
    ds_out = xr.Dataset({
        # 'time': ( ['time'], dsA['time'] ),
        'lat': (['lat'], grids[grid2use]['lat']),
        'lon': (['lon'], grids[grid2use]['lon']),
    })
    # Create a regidder (to be reused )
    regridder = xe.Regridder(dsA, ds_out, method, reuse_weights=True)
    # loop and regrid variables
    ds_l = []
    vars2regrid = dsA.data_vars
    for var2use in vars2regrid:
        # create a dataset to re-grid into
        ds_out = xr.Dataset({
            'lat': (['lat'], grids[grid2use]['lat']),
            'lon': (['lon'], grids[grid2use]['lon']),
        })
        # get a DataArray
        dr = dsA[var2use]
        # build regridder
        dr_out = regridder(dr)
        # Important note: Extra dimensions must be on the left, i.e. (time, lev, lat, lon) is correct but (lat, lon, time, lev) would not work. Most data sets should have (lat, lon) on the right (being the fastest changing dimension in the memory). If not, use DataArray.transpose or numpy.transpose to preprocess the data.
        # exactly the same as input
#        xr.testing.assert_identical(dr_out['time'], dsA['time'])
        # save variable
        ds_l += [dr_out]
    ds = xr.Dataset()
    for n, var2use in enumerate(vars2regrid):
        ds[var2use] = ds_l[n]
    # clean up and return dataset
    regridder.clean_weight_file()
    return ds


def add_raster_of_country2ds(ds, country='South Africa',
                             set_all_regions2one=True,
                             test_plot=False, dpi=320):
    """
    Add raster outline of country to spatial dataset
    """
    # Get shapes for country
    shapes = get_shapes4country(country=country)
    # Add country's states as a layer
    ds['states'] = rasterize(shapes, ds.coords)
    # Test plot of this?
    if test_plot:
        from . plotting import quick_map_plot
        savename = 'spatial_plot_of_shapes4country_{}'.format(country)
        quick_map_plot(ds, var2plot='states', savename=savename)

    # set all the regions (e.g. counties/states) in a country to 1
    if set_all_regions2one:
        arr = ds['states'].values
        arr[np.where(~np.isnan(arr))] = 1
        ds['states'].values = arr
    return ds


def get_shapes4country(country='South Africa'):
    """
    Get shapes (polygons) for country from Natural earth

    NOTES
     - data credit    http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
    """
    # location of data
    URL = "http://www.naturalearthdata.com/downloads/10m-cultural-vectors"
    URL += "/10m-admin-1-states-provinces/"
    # Shapefiles locally?
    # TODO - update to download automatically and store in AC_tools' data directory
    shapefiles = 'ne_10m_admin_1_states_provinces_lakes'
#    shapefiles = 'ne_10m_admin_1_states_provinces'
    folder = '/mnt/lustre/users/ts551/labbook/Python_progs/'
    folder += '/AC_tools/data/shapefiles/{}'.format(shapefiles, shapefiles)
    states = geopandas.read_file(folder)
    # Just select state of interest
    choosen_states = states.query("admin == '{}'".format(country))
    choosen_states = choosen_states.reset_index(drop=True)
    # Get the shapes
    shapes = zip(choosen_states.geometry, range(len(choosen_states)))
    return shapes


def transform_from_latlon(lat, lon):
    """
    Tranform from latitude and longitude
    NOTES:
     - credit - Shoyer https://gist.github.com/shoyer/0eb96fa8ab683ef078eb
    """
    from affine import Affine
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale


def rasterize(shapes, coords, fill=np.nan, **kwargs):
    """
    Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.

    NOTES:
     - credit - Shoyer https://gist.github.com/shoyer/0eb96fa8ab683ef078eb
    """
    from rasterio import features
    transform = transform_from_latlon(coords['lat'], coords['lon'])
    out_shape = (len(coords['lat']), len(coords['lon']))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    return xr.DataArray(raster, coords=coords, dims=('lat', 'lon'))


def ocean_unmasked(res='4x5', debug=False):
    """
    Get ocean mask from GEOS-Chem LWI ( at given resolution)

    NOTES:
     - this currently returns mask with all except ocean masked
     - NEEDS UPDATE
    """

    from .GEOSChem_bpch import get_LWI_map
    if debug:
        print(('ocean_mask called for: ', res))

    # Create a mask from land/water/ice indices
    m = np.ma.masked_not_equal(get_LWI_map(res=res), 0)
    if debug:
        print((mask, mask.shape))
    return m.mask


def land_unmasked(res='4x5', debug=False):
    """
    Get land mask from GEOS-Chem LWI ( at given resolution)
    """
    from .GEOSChem_bpch import get_LWI_map  # Kludge, use GEOS-Chem LWI

    # Create a np.ma mask
    if debug:
        print(('land_mask called for: ', res))
    m = np.ma.masked_not_equal(get_LWI_map(res=res), 1)
    if debug:
        print((mask, mask.shape))
    return m.mask


def ice_unmasked(res='4x5', debug=False):
    """
    Get mask for all areas but ice surface. Ice mask is from GEOS-Chem
    LWI ( at given resolution)

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask
    m = np.logical_not((land_unmasked(res)*ocean_unmasked(res)))
    if debug:
        print((mask, mask.shape))
    return m


def surface_unmasked(res='4x5', trop_limit=False, mask2D=False,
                     debug=False):
    """
    Get mask for all areas but surface

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask
    m = np.ma.array(np.zeros(get_dims4res(res)), mask=False)
    m[..., 0] = 1
    m = np.ma.masked_not_equal(m, 1)

    if trop_limit:
        m = m[..., :38]
    if debug:
        print((mask, mask.shape))

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def tropics_unmasked(res='4x5', saizlopez=False, pradosroman=False,
                     mask2D=False):
    """
    Get mask for all areas but tropics

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    if saizlopez or pradosroman:
        lats = np.arange(-20, 20, 1)
    else:
        lats = np.arange(-22, 22, 1)
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def mid_lats_unmasked(res='4x5', saizlopez=False, pradosroman=False,
                      mask2D=False):
    """
    Get mask for all areas but mid latitudes

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    if saizlopez:
        lats = np.concatenate((np.arange(-50, -20, 1), np.arange(20, 50, 1)))
    elif pradosroman:
        lats = np.concatenate((np.arange(-60, -20, 1), np.arange(20, 61, 1)))
    if ((not saizlopez) and (not saizlopez) and (res == '4x5')):
        lats = np.concatenate((np.arange(-50, -26, 1), np.arange(26, 51, 1)))
    if ((not saizlopez) and (not saizlopez) and (res == '2x2.5')):
        lats = np.concatenate((np.arange(-50, -24, 1), np.arange(24, 51, 1)))

    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def mask_lat40_2_40(res='4x5', mask2D=False):
    """
    Get mask for all areas but 40S-40N latitudes (40 to 40 deg latitude)

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    lats = np.arange(-42, 42, 1)  # use 42 due to 4x5 grid
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def extratropics_unmasked(res='4x5', mask2D=False):
    """ Get mask for all areas but extratropics
    NOTES:
        - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    lats = np.concatenate((np.arange(-89, -26, 1), np.arange(26, 90, 1)))
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1
    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def all_unmasked(res='4x5', mask2D=False):
    """ Get un-masked mask of size GEOS-Chem dimensions
    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    m = np.ma.array(np.ones(get_dims4res(res)), mask=False)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def mask_3D(hPa, sect, MBL=True, res='4x5', extra_mask=None,
            M_all=False, use_multiply_method=True, trop_limit=False,
            verbose=True, debug=False):
    """
    Creates Maskes by pressure array  (required shape: 72,46,47),
    with conditions (lower and upper bounds) set by given cases for
    MBL, UT, FT

    Parameters
    -------
    sect (Str): section of the atmosphere of interest (e.g. MBL, UT...)
    hPa (array): array for pressures ( in hPa)
    MBL (bool): apply a mask for the marine boundary layer
    res (str): the resolution of required output/input arrays (e.g. '4x5' )
    use_multiply_method (bool): Create arrays of ones and zeros
    trop_limit (bool): limit 3D arrays to troposphere
    debug (bool): legacy debug option, replaced by python logging
    verbose (bool): legacy debug option, replaced by python logging
    extra_mask (str): name of additional region (e.g. ocean) to mask
    M_all (bool): apply oceanic masking to all regions

    Returns
    -------
    (np.ma.mask)

    NOTES:
     - originally written to generate masks for mulitplication
    (i.e. use_multiply_method = True ), but can also be use to make
    more pythonic masks ( if use_multiply_method=False )
    """
    if verbose:
        print(('mask_3D called for sect={}, use_multiply_method={}'.format(
            sect, use_multiply_method) + ', M_all={}, and '.format(M_all) +
            'with debug: {}, verbose:{}'.format(sect, debug, verbose)))

    # Get atmospheric region as case defining lower and upper bounds
    cases = {
        'BL': [1200., 900.], 'MBL': [1200., 900.], 'FT': [900., 350.],
        'UT': [350., 75.], 'All': [1200., 75.]
    }
    l, h = cases[sect]

    # --- Mask between upper and lower values
    m = np.ones(get_dims4res(res))
    m[(hPa >= l)] = 0
    m[(hPa < h)] = 0
    logging.debug('Sect={}, l={}, h={}'.format(sect, l, h))
    logging.debug('{}'.format(
        *[[i.min(), i.max(), i.mean(), i.sum(), i.shape] for i in [m]]))

    # Mask off the 'sect' area that still equals 1
    m = np.ma.masked_equal(m, 1)

    # --- Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        m = m[..., :38]

    if not isinstance(extra_mask, type(None)):
        # only consider MBL
        if (MBL and sect == 'BL') or (sect == 'MBL') or M_all:
            if use_multiply_method:
                return m.mask * extra_mask * land_unmasked(res)
            else:
                print("WARNING: needs 3D arrays, use 'mask_all_but' instead")
                sys.exit()
        else:
            if use_multiply_method:
                return m.mask * extra_mask
            else:
                print("WARNING: needs 3D arrays, use 'mask_all_but' instead")
                sys.exit()

    # --- Only consider MBL (or MFT/MFT for Saiz-Lopez 2014 comparison)
    if (MBL and sect == 'BL') or (sect == 'MBL') or M_all:
        if use_multiply_method:
            return m.mask * land_unmasked(res)
        else:
            land_unmasked_ = mask_all_but('Land', mask3D=True, res=res,
                                          use_multiply_method=False,
                                          trop_limit=trop_limit)

            # MBL unmasked
            m = np.logical_or(np.logical_not(m.mask),
                              np.logical_not(land_unmasked_))

            # Invert as function expected to return oposite
            m = np.logical_not(m)

            return m   # m is a mask here

    return m.mask


def lat2lat_2D_unmasked(lowerlat=None, higherlat=None, res='2x2.5',
                        debug=False):
    """
    Takes a lower and higher latitude value and then creates mask of area boxes outside
    given laitidue limits.
    """

    # Get vars
    lon_c, lat_c, NIC = get_latlonalt4res(res=res, centre=True)

    # mask between upper and lower values
    lats = [i for i in lat_c if ((i >= lowerlat) and (i < higherlat))]
    lats = [get_gc_lat(i, res=res) for i in lats]

    # fill all lat and lon True or False
    m = np.zeros(get_dims4res(res))[:, :, 0]
    print((m.shape, np.sum(m)))
    for i in lats:
        m[:, i] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask


def southpole_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with southpole masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5':
        lats = np.arange(-89, -62, 1)  # define S pole as > 60S
    else:
        lats = np.arange(-89, -60, 1)  # define S pole as > 60S
#    lats = np.arange(-89, -80,1 ) # define S pole as > 80S
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def northpole_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with northpole masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of zeros
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5':
        lats = np.arange(62, 90, 1)  # define N pole as > 60N
    else:
        lats = np.arange(60, 90, 1)  # define N pole as > 60N
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def NH_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with north hemisphere masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    if res == '4x5':
        lats = np.arange(1, 91, 1)
    elif res == '2x2.5':
        lats = np.arange(0, 89, 1)
        print('CHECK (NH) mask for non 4x5 resolutions')
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :].mask = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def SH_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with south hemisphere masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    if res == '4x5':
        lats = np.arange(-89, 0, 1)
    if res == '2x2.5':
        lats = np.arange(-90, 0, 1)
        print('CHECK (SH) mask for non 4x5 resolutions')
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :].mask = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def location_unmasked(res='4x5', lat=None, lon=None, mask2D=False):
    """
    Return global np.ma.mask with all locations apart from that containing the
    location (lat, lon) masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    # Get location index
    assert all([type(i) == float for i in (lat, lon)]
               ), 'lat & lon must be floats'
    lat_ind = get_gc_lat(lat=lat, res=res)
    lon_ind = get_gc_lon(lon=lon, res=res)
    # Unmask grid box for location
    m.mask[lon_ind, lat_ind] = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


def get_analysis_masks(masks='basic',  hPa=None, M_all=False, res='4x5',
                       saizlopez=False, r_pstr=True, wd=None, trop_limit=True,
                       mask4D=False,
                       use_multiply_method=True, debug=False):
    """
    Return list of mask arrays for analysis

    NOTES:
    - For comparisons with Saiz-Lopez et al. 2014, set M_all to True,
    this masks the UT and FT as well as the BL (therefore gives MBL, MFT, MUT)

    """
    # --- Check hPa has been provided as arg.
    if isinstance(hPa, type(None)):
        print('ERROR: Please provide array of hPa to get_analysis_masks')

    if masks == 'full':
        # ---- List of masks
        mtitles = [
            'All', 'Ocean', 'Land', 'Ice', 'All Sur.',  'Ocean Sur.',
            'Land Sur.', 'Ice Sur.', 'NH', 'SH', 'Tropics', 'Ex. Tropics',
            'Mid Lats', 'Ocn. 50S-50N',  '50S-50N'
        ]
        tsects3D = ['MBL', 'BL', 'FT',  'UT']
        # get MBL, FT and UT maskes
        sects3D = [
            mask_3D(hPa, i, MBL=False, M_all=M_all, res=res,
                    use_multiply_method=use_multiply_method) for i in tsects3D
        ]
        # ---- Use non-pythonic mulitply method?
        if use_multiply_method:
            maskes = [mask_all_but(i, trop_limit=trop_limit, mask3D=True,
                                   use_multiply_method=True, res=res)
                      for i in mtitles]
            # if comparison with saiz-lopez 2014,
            if M_all:
                ind = [n for n, i in enumerate(mtitles) if not ('MBL' in i)]
                for n in ind:
                    maskes[n] = maskes[n]*land_unmasked(res=res)
        # --- Use pythonic approach
        else:
            maskes = [mask_all_but(i, trop_limit=trop_limit, mask3D=True,
                                   use_multiply_method=False, res=res)
                      for i in mtitles]
            # If not 'use_multiply_method', then invert hPa masks
            sects3D = [np.logical_not(i) for i in sects3D]
        # Add to mask and mask title lists
        maskes = maskes + sects3D
        mtitles = mtitles + tsects3D
        # Also create print strings...
        npstr = '{:<12}'*len(maskes)
        pstr = '{:<12,.2f}'*len(maskes)
    if masks == 'basic':
        tsects3D = ['All', 'MBL', 'BL', 'FT',  'UT']
        mtitles = [i+' (All)' for i in tsects3D] + \
            [i+' (Tropics)' for i in tsects3D] +   \
            [i+' (Mid Lats)' for i in tsects3D]
        # Standard maskes none, tropics, mid-lats (3)
        maskes = [
            np.logical_not(i) for i in (all_unmasked(res=res),
                                        tropics_unmasked(
                                            res, saizlopez=saizlopez),
                                        mid_lats_unmasked(res))]
        # Additional masks - tsects3D (4+1) * standard maskes (3)
        dmaskes = [
            [mask_3D(hPa, i, MBL=False, extra_mask=mask, M_all=M_all, res=res)
             for i in tsects3D] for mask in maskes]
        # Unpack and set as maskes (12) to list (3)
        dmaskes = [j for i in dmaskes for j in i]
        print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        maskes = dmaskes
        print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        # If comparison with saiz-lopez 2014 appli marine mask to all...
        if M_all:
            ind = [n for n, i in enumerate(mtitles) if not 'MBL' in i]
            for n in ind:
                maskes[n] = maskes[n]*land_unmasked(res=res)
        if debug:
            print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        # Also create print strings...
        npstr = '{:<15}'*len(maskes)
        pstr = '{:<15,.2f}'*len(maskes)
    if masks == 'trop_regions':
        mtitles = ['BL', 'FT',  'UT']
        maskes = [mask_3D(hPa, i, M_all=M_all, MBL=False, res=res)[:, :, :38]
                  for i in mtitles]
        # Also create print strings...
        npstr = '{:<15}'*len(maskes)
        pstr = '{:<15,.2f}'*len(maskes)
    # Only consider the "chemical troposphere" - according v9-2
    if trop_limit:
        maskes = [i[:, :, :38] for i in maskes]
    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        for n, mask in enumerate(maskes):
            if any([(mask.shape[-1] == i) for i in [12]]):
                pass
            else:  # concatenate dimensions
                maskes[n] = np.concatenate([mask[..., None]]*12, axis=3)
    if r_pstr:
        return maskes, mtitles, npstr, pstr
    else:
        return maskes, mtitles


def mask_all_but(region='All', M_all=False, saizlopez=False,
                 res='4x5', trop_limit=True, mask2D=False, mask3D=False,
                 mask4D=False,
                 use_multiply_method=True, lat=None, lon=None,
                 verbose=False, debug=False):
    """
    Mask selector for analysis. global mask provided for with given region
        unmasked

    Parameters
    -------
    res (str): the resolution if wd not given (e.g. '4x5' )
    M_all (bool): maask all marine areas?
    saizlopez (bool): use tropics definition from Saiz-Lopez er al 2014
    trop_limit (bool): limit 4D arrays to troposphere
    mask2D/mask3D/mask4D(booolean): ensure mask returned is 2D/3D/4D
    use_multiply_method (bool): return array of ones, that can be mulitpled
    through an array to set data to zero
    verbose (bool): legacy debug option, replaced by python logging
    debug (bool): legacy debug option, replaced by python logging
    loc (str): location
    lon, lat (float): lat/lon locations to leave nearest grid box unmasked

    Returns
    -------
    (np.ma.mask) or (np.array) (later if use_multiply_method==True)

    Notes
    -----
    "unmask_all" yeilds completely unmasked array
    function was oringialyl used to mulitple masks, however, this approch is
    unpythonic and therefore reccomended against.
    """
    logging.info('mask_all_but called for region {}'.format(region))
    # --- Setup cases...
    # ( except None, unmask_all and global to retrive no mask )
    case = {
        'Tropics': 0,
        'tropics': 0,
        'mid_lats': 1,
        'Mid Lats': 1,
        'Mid lats': 1,
        'south_pole': 2,
        'south pole': 2,
        'north_pole': 3,
        'north pole': 3,
        None: 4,
        'unmask_all': 4,
        'All': 4,
        'global': 4,
        # NEED TESTING ...
        'Extratropics': 5,
        'Ex. Tropics': 5,
        'Oceanic': 6,
        'Ocean': 6,
        'NH': 7,
        'SH': 8,
        'Ice': 10,
        'Land': 11,
        'lat40_2_40': 12,
        'Ocean Tropics': 13,
        'Oceanic Tropics': 13,
        'Ocn. Trop.': 13,
        'Land Tropics': 14,
        'All Sur.': 15,
        'surface': 15,
        'Ocean Sur.': 16,
        'Land Sur.': 17,
        'Ice Sur.': 18,
        'lat50_2_50': 19,
        '50S-50N': 19,
        #    'Oceanic lat50_2_50': 20,
        'Ocn. 50S-50N': 20,
        #     'South >60': 2,
        #      'North >60': 3
        'North Sea': 21,
        'Med. Sea': 22,
        'Mediterranean Sea': 22,
        'Black Sea': 23,
        'Irish Sea': 24,
        'Europe': 25,
        'EU': 25,
        #    'Surface BL': 26,
        'Land Tropics Sur.': 27,
        'Boreal Land': 28,
        'Alps':  29,
        'loc': 30,
        'location': 30,
        'France': 31,
    }[region]

    # --- This is a simple way of using masks ( as multiplers )
    # i.e. all (future) functions should have use_multiply_method=False
    # and not use the code below
    if use_multiply_method:  # Kludge
        print(('!'*50, 'WARNING: using mulitply method for masking. '))

        # For case, pull mask from case list
        if case == 0:
            mask = tropics_unmasked(res=res, saizlopez=saizlopez)
        elif case == 1:
            mask = mid_lats_unmasked(res=res)
        elif case == 2:
            mask = southpole_unmasked(res=res)
        elif case == 3:
            mask = northpole_unmasked(res=res)
        elif case == 4:
            #        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked(res=res)
        elif case == 5:
            mask = extratropics_unmasked(res=res)
        elif case == 6:
            mask = ocean_unmasked(res=res)
        elif case == 7:
            mask = NH_unmasked(res=res)
        elif case == 8:
            mask = SH_unmasked(res=res)
        elif case == 10:
            mask = ice_unmasked(res=res)
        elif case == 11:
            mask = land_unmasked(res=res)
        elif case == 12:
            mask = mask_lat40_2_40(res=res)
        elif case == 13:  # 'Oceanic Tropics'
            mask = np.ma.mask_or(ocean_unmasked(res=res),
                                 tropics_unmasked(res=res,
                                                  saizlopez=saizlopez))
        elif case == 14:  # 'Land Tropics'
            mask = np.ma.mask_or(land_unmasked(res=res),
                                 tropics_unmasked(res=res,
                                                  saizlopez=saizlopez))
        elif case == 15:  # 'All Sur.'
            mask = surface_unmasked(res=res)
        elif case == 16:  # 'Ocean Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ocean_unmasked(res=res))
        elif case == 17:  # 'Land Sur.':
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 land_unmasked(res=res))
        elif case == 18:  # 'Ice Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ice_unmasked(res=res))
        elif case == 19:  # '50S-50N'
            mask = lat2lat_2D_unmasked(lowerlat=-50, higherlat=50,
                                       res=res)
        elif case == 20:  # 'Ocn. 50S-50N'
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=-50,
                                                     higherlat=50, res=res),
                                 ocean_unmasked(res=res)[..., 0])
        elif case == 21:
            mask = get_north_sea_unmasked(res=res)
        elif case == 25:
            mask = get_EU_unmasked(res=res)
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )
        elif case == 27:  # 'Land Tropics Sur.':
            tmp = np.ma.mask_or(surface_unmasked(res=res),
                                land_unmasked(res=res))
            mask = np.ma.mask_or(tmp, tropics_unmasked(res=res))
        else:
            print('WARNING - Mask not setup for case={}'.format(case))
            sys.exit()
        # Invert mask to leave exception unmasked if used to multiply
        mask = np.logical_not(mask)

    # --- This is a more pythonic way of using masks (Use as preference)
    else:
        # For case, pull mask from case list
        if case == 0:
            mask = tropics_unmasked(res=res, saizlopez=saizlopez)
        elif case == 1:
            mask = mid_lats_unmasked(res=res)
        elif case == 2:
            mask = southpole_unmasked(res=res)
        elif case == 3:
            mask = northpole_unmasked(res=res)
        elif case == 4:
            #        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked(res=res)
        elif case == 5:
            mask = extratropics_unmasked(res=res)
        elif case == 6:
            mask = ocean_unmasked(res=res)
        elif case == 7:
            mask = NH_unmasked(res=res)
        elif case == 8:
            mask = SH_unmasked(res=res)
        elif case == 10:
            mask = ice_unmasked(res=res)
        elif case == 11:
            mask = land_unmasked(res=res)
        elif case == 12:
            mask = mask_lat40_2_40(res=res)
        elif case == 13:
            mask = np.ma.mask_or(ocean_unmasked(res=res),
                                 tropics_unmasked(res=res,
                                                  saizlopez=saizlopez))
        elif case == 14:
            mask = np.ma.mask_or(land_unmasked(res=res),
                                 tropics_unmasked(res=res,
                                                  saizlopez=saizlopez))
        elif case == 15:  # 'All Sur.'
            mask = surface_unmasked(res=res)
        elif case == 16:  # 'Ocean Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ocean_unmasked(res=res))
        elif case == 17:  # 'Land Sur.':
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 land_unmasked(res=res))
        elif case == 18:  # 'Ice Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ice_unmasked(res=res))
        elif case == 19:
            mask = lat2lat_2D_unmasked(lowerlat=-50, higherlat=50,
                                       res=res)
        elif case == 20:
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=-50,
                                                     higherlat=50, res=res),
                                 ocean_unmasked(res=res)[..., 0])
        elif case == 21:
            mask = get_north_sea_unmasked(res=res)
        elif case == 22:
            mask = get_mediterranean_sea_unmasked(res=res)
        elif case == 23:
            mask = get_unmasked_black_sea(res=res)
        elif case == 24:
            mask = get_unmasked_irish_sea(res=res)
        elif case == 25:
            mask = get_EU_unmasked(res=res)
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )
        elif case == 27:  # 'Land Tropics Sur.':
            tmp = np.ma.mask_or(surface_unmasked(res=res),
                                land_unmasked(res=res))
            mask = np.ma.mask_or(tmp, tropics_unmasked(res=res))
        elif case == 28:
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=50,
                                                     higherlat=80, res=res),
                                 land_unmasked(res=res)[..., 0])
        elif case == 29:  # Alps
            # Alps mask
            lowerlat = 43
            higherlat = 47
            lowerlon = 5
            higherlon = 15
            # Get a mask for lat and lon range, then combine
            mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                        higherlat=higherlat)
            mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                        higherlon=higherlon)
            mask = np.ma.mask_or(mask1, mask2)
        elif case == 30:  # Location ('loc' )
            # Alps
            mask = location_unmasked(lat=lat, lon=lon, res=res)
        elif case == 31:  # Rough(!) France map
            # mask
            mask = get_France_unmasked(res=res)

        else:
            print('WARNING - Mask not setup for case={}'.format(case))
            sys.exit()

    logging.debug('prior to setting dimensions: {}'.format(mask.shape))
    # Apply Saiz-Lopez Marine MFT/MUT? <= should this be before multiply op.?
    if M_all:
        if use_multiply_method:  # Kludge
            mask = mask*land_unmasked(res=res)
        else:
            # check this!!!
            mask = np.ma.mask_or(mask, land_unmasked(res=res))

    # Ensure returned arrays are 2D
    if mask2D:
        if len(mask.shape) == 2:
            pass
        elif len(mask.shape) == 3:
            mask = mask[..., 0]
        elif len(mask.shape) == 4:
            mask = mask[..., 0, 0]

    # Create 3D array by concatenating through altitude dimension
    if mask3D:
        if any([(mask.shape[-1] == i) for i in (38, 47)]):
            pass
        else:  # concatenate dimensions
            if len(mask.shape) == 3:
                mask = np.concatenate([mask]*47, axis=2)
            elif len(mask.shape) == 2:
                mask = np.concatenate([mask[..., None]]*47, axis=2)

    # Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        if (len(mask.shape) == 2) or mask2D:
            pass
        else:
            mask = mask[..., :38]

    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        if any([(mask.shape[-1] == i) for i in [12]]):
            pass
        else:  # concatenate dimensions
            mask = np.concatenate([mask[..., None]]*12, axis=3)
    logging.debug('post to setting dimensions: {}'.format(mask.shape))
    logging.info("returning a 'mask' of type:{}".format(type(mask)))
    return mask


def lon2lon_2D_unmasked(lowerlon, higherlon, res='2x2.5', debug=False):
    """
    Takes a lower and higher latitude value and then creates
    mask to given given limits.
    """

    # Get vars
    lon_c, lat_c, NIU = get_latlonalt4res(res=res, centre=True)
    if debug:
        print((lon_c, lowerlon, higherlon))

    # Mask between upper and lower values
    lons = [i for i in lon_c if ((i >= lowerlon) and (i < higherlon))]
    lons = [get_gc_lon(i, res=res) for i in lons]

    # Fill all lat and lon True or False
    m = np.zeros(get_dims4res(res))[:, :, 0]
    if debug:
        print((m.shape, np.sum(m)))
    for i in lons:
        m[i, :] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask


def get_EU_unmasked(res='1x1'):
    """
    Mask 'EU' as defined by GEOS-Chem EU grid the grid of "'0.5x0.666" resolution is
    used by default, but any list of lat and lons could be provided and the extremities
    would be used as the mask edges
    """
    EU_resolutions = ['0.25x0.3125', '0.5x0.666']
    if res not in EU_resolutions:
        EU_res = EU_resolutions[0]  # default = '0.25x0.3125'
#        EU_res = EU_resolutions[1] # default = '0.5x0.666'
    else:
        EU_res = res
    # Get GEOS-Chem EU lat and lons
    lon, lat, NIU = get_latlonalt4res(res=EU_res)
    # mask lats
    m1 = lat2lat_2D_unmasked(lowerlat=lat.min(), higherlat=lat.max(), res=res)
    # Mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=lon.min(), higherlon=lon.max(), res=res)
    # Combine maskes
    m = m1 + m2
    return m


def get_cruise_track_mask(max_lon=None, min_lon=None, max_lat=None,
                          min_lat=None, unmask_water=True, res='4x5',
                          trop_limit=True):
    """
    Mask whole area of ship based research campaigns for bulk comparison
    """
    # only look at surface
    m = surface_unmasked(res=res, trop_limit=trop_limit)
    # apply ocean mask
    if unmask_water:
        m = m + ocean_unmasked(res=res)
    # Mask over given longitude range, if provided
    if not isinstance(max_lon, type(None)):
        m = m + lon2lon_2D_unmasked(lowerlon=min_lon, higherlon=max_lon,
                                    res=res)[:, :, None]
    # Mask over given latitude range, if provided
    if not isinstance(max_lat, type(None)):
        m = m + lat2lat_2D_unmasked(lowerlat=min_lat, higherlat=max_lat,
                                    res=res)[:, :, None]
    # Invert
    m = np.logical_not(m)
    return m


def get_France_unmasked(res='4x5'):
    """ A rough Mask of France for use with 2x2.5 / 4x5 model output """
    # France mask
    lowerlat = 42.5
    higherlat = 51
    lowerlon = -4.441
    higherlon = 7.7577
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Only consider land grid boxes
    mask = np.ma.mask_or(mask, land_unmasked(res=res)[..., 0])
    return mask


def get_mediterranean_sea_unmasked(res='0.25x0.3125'):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 34
    lowerlon = -6.5
    # East North corner (~Ukraine)
    higherlat = 47
    higherlon = 38
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    # Also remove black sea ( by removing an inverted unmasked mask )
    mask3 = get_black_sea_unmasked(res=res, unmask_water=False)
    mask = np.ma.mask_or(mask, np.logical_not(mask3))
    # Also remove bay of biscay
    # Also remove black sea
    mask4 = get_bay_of_biscay_unmasked(res=res)
    mask = np.ma.mask_or(mask, np.logical_not(mask4))
    return mask


def get_black_sea_unmasked(res='0.25x0.3125', unmask_water=True):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 41
    lowerlon = 26.8
    # East North corner (~Ukraine)
    higherlat = 50
    higherlon = 43
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    if unmask_water:
        mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    return mask


def get_bay_of_biscay_unmasked(res='0.25x0.3125', unmask_water=True):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 42.5
    lowerlon = -10
    # East North corner (~Ukraine)
    higherlat = 51
    higherlon = 0
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    if unmask_water:
        mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    # Also remove black sea
    return mask


def get_Southern_Africa_Masked(res='0.25x0.3125', ):
    """
    Create a mask of sub equatorial africa
    """
    # Sub equatorial africa
    lowerlat = -40
    higherlat = 0
    lowerlon = -10
    higherlon = 55
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)

    return mask


def get_CVAO_Africa_nest_Masked(res='0.25x0.3125', ):
    """
    Create a mask of sub equatorial africa
    """
    # Sub equatorial africa
    lowerlat = 0.0
    higherlat = 34.0
    lowerlon = -32.0
    higherlon = 15.0
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)

    return mask


# def get_north_sea_unmasked( res='0.25x0.3125' ):
#     """
#     A rough Mask of the North Sea for use with ~0.5/~0.25 mdodel output.
#     (inc. English channel. )
#
#     """
#     mask latitudes of north sea
#     Drawing a box that includes all of North sea and English channel
#     Near Brest in France 48.3669927,-4.7560745
#     Lillehammer 61.1122408,10.4386779
#
#     mask lats
#     m1 = lat2lat_2D_unmasked( lowerlat=48.3669927, higherlat=61.1122408, res=res )
#
#     mask lons
#     m2 = lon2lon_2D_unmasked(lowerlon=-4.7560745, higherlon=10.4386779, res=res )
#
#      combine maskes
#     m = m1 + m2
#
#     remove all water.
#    m = np.ma.mask_or( m, ocean_unmasked( res=res)[...,0] ) # old approach
#     the below will work for muliple options.
#    m = m +
#    m = ocean_unmasked( res=res)[...,0]
#
#
#     remove irish sea
#
#     return m.mask
#
#
# def get_unmasked_black_sea( res='0.25x0.3125'):
#     """
#     A rough Mask of the Black Sea for use with ~0.5/~0.25 mdodel output.
#     """
#
#     Drawing a box that include all of Black. Sea
#     East South corner (south Jordan) = 43.4623268,33.3809392
#     West North corner  = 44.17207,25.30604
#     pass
#
#
# def get_unmasked_irish_sea( res='0.25x0.3125', unmasked_oceans=None):
#     """
#     A rough Mask of the Irish Sea for use with ~0.5/~0.25 mdodel output.
#     """
#
#     NE corner Glasgow - 55.8553803,-4.3725463
#     SW corner Cork - 51.8959842,-8.5332609
#     mask lats
#     m1 = lat2lat_2D_unmasked( lowerlat=51.8959842, higherlat=55.855380, res=res )
#
#     mask lons
#     m2 = lon2lon_2D_unmasked(lowerlon=8.5332609, higherlon=4.3725463, res=res )
#
#      combine maskes
#     m = m1 + m2
#
#     only consider oceans
#     m = np.ma.mask_or( ocean_unmasked( res=res)[...,0], m )
#
#
#     return mask

def get_2D_nighttime_mask4date_pd(date=None, ncfile=None, res='4x5',
                                  mask_daytime=False, buffer_hours=0,
                                  debug=False):
    """
    Creates 2D (lon,lat) masked (1=Masked) for nighttime for a given list of
    dates

    Parameters
    -------
    date (datetime): date to use (UTC)
    mask_daytime (bool): mask daytime instead of nightime
    ncfile (str): location to netCDF file - not implemented...
    res (str): resolution, if using resolutions listed in get_latlonalt4res
    buffer_hours (float/int): number of hours to buffer subrise/sunset with
     (This will act to increase the size of the mask - e.g. if masking
      nightime, then an extra hour of nightime would be added to sunrise, and
      removed from sunset. )

    Returns
    -------
    (np.array) with conditional values masked (1=masked)

    ncfile (NetCDF file): NetCDF file to extract lat and lon metadata from

    Notes
    -----
     - if ncfile provide programme will work for that grid.
    """
    # Astronomical math
    import ephem
    from ephem import AlwaysUpError, NeverUpError
    # And functions in other AC_tools modules
    from .AC_time import add_days, add_hrs
    logging.info('get_2D_nighttime_mask4date_pd called for {}'.format(date))

    # Profile function...
    if debug:
        start_time = time.time()

    # --- Local variables?
    # reference data for ephem (number of days since noon on 1899 December 31)
    ref_date = datetime.datetime(1899, 12, 31, 12)

    # --- Get LON and LAT variables
    if isinstance(ncfile, type(None)):
        # extract from refence files
        lons, lats, alts = get_latlonalt4res(res=res)
    else:
        # TODO - allow any lat, lon grid to be used by taking input lats and
        # lons from ncfile file/arguments.
        print('Not implemented')
        sys.exit()
    if debug:
        print(("--- (start-1) %s seconds ---" % (time.time() - start_time)))

    # --- setup function to mask based on date, lat and lon
    def mask_nighttime(lon, lat, date=date, mask_daytime=mask_daytime,
                       ref_date=datetime.datetime(1899, 12, 31, 12),
                       buffer_hours=buffer_hours, debug=False):
        """
        sub-function to mask if nightime for a given date at a specific lat/lon
        """
        # --- get lat and lon values from columns
        if debug:
            print(("--- (s4-1) %s seconds ---" % (time.time() - start_time)))
        # --- get sunrise and sunset for location
        o = ephem.Observer()
        # set lat (decimal?), lon (decimal?), and date (UTC)
        o.lat = str(lat)
        o.long = str(lon)
        o.date = date
        # planetary body
        s = ephem.Sun()
        if debug:
            print(("--- (s4-2) %s seconds ---" % (time.time() - start_time)))

        # Compute sun vs observer
        s.compute()
        if debug:
            print(("--- (s4-3) %s seconds ---" % (time.time() - start_time)))

        # Work out if day or night based on sunrises and sunsets
        mask_value = 0
        try:

            # get sunrise time and date
            next_rising = o.next_rising(s)
            next_setting = o.next_setting(s)

            # convert to datetime.datetime
            next_rising = add_days(ref_date, next_rising)
            next_setting = add_days(ref_date, next_setting)

            # did the sun last rise or set? (inc. any adjustments)
            sun_last_rose = False
            if next_setting < next_rising:
                sun_last_rose = True

            # Add buffer to rising/setting if provided with buffer_hours
            if buffer_hours != 0:

                # Calculate last rise
                previous_rising = o.previous_rising(s)
                # convert to datetime.datetime
                previous_rising = add_days(ref_date, previous_rising)
                # Calculate last setting
                previous_setting = o.previous_setting(s)
                # convert to datetime.datetime
                previous_setting = add_days(ref_date, previous_setting)

                # Calculate absolute difference
                time_from_rise = (date-previous_rising).total_seconds()
                time_till_set = (date-next_setting).total_seconds()
                time_from_set = (date-previous_setting).total_seconds()
                time_till_rise = (date-next_rising).total_seconds()

                # If absolutely difference less than buffer
                if abs(time_from_rise)/60./60. <= buffer_hours:
                    mask_value = 1
                elif abs(time_till_set)/60./60. <= buffer_hours:
                    mask_value = 1
                elif abs(time_from_set)/60./60. < buffer_hours:
                    mask_value = 1
                elif abs(time_till_rise)/60./60. < buffer_hours:
                    mask_value = 1

            # --- Check if daytime or nighttime and mask if condition met.
            if sun_last_rose:
                if mask_daytime:
                    # ... and has not set yet, it must be daytime
                    if (date < next_setting):
                        mask_value = 1

            # if the sun last set... (mask nighttime is default)
            else:
                # if mask nighttime (aka not mask_daytime)
                if not mask_daytime:
                    # ... and has not risen yet, it must be nighttime
                    if (date < next_rising):
                        mask_value = 1

        # Add gotcha for locations where sun is always up.
        except AlwaysUpError:
            if mask_daytime:
                mask_value = 1

        # Add gotcha for locations where sun is always down.
        except NeverUpError:
            if not mask_daytime:
                mask_value = 1

        except:
            print('FAIL')
            sys.exit()

        # Mask value in array
        return mask_value

    # --- Setup an unstack(ed) pandas dataframe to contain masked values
    if debug:
        print(("--- (2) %s seconds ---" % (time.time() - start_time)))
    # Use list comprehension to setup list of indices for lat and lon
    # Better way of doing this? (e.g. pd.melt?)
    ind_lat_lons_list = [[lon_, lat_] for lat_ in lats for lon_ in lons]
    if debug:
        print(("--- (3) %s seconds ---" % (time.time() - start_time)))
    # Make this into a pd.DataFrame and label columns.
    df = pd.DataFrame(ind_lat_lons_list)
    df.columns = ['lons', 'lats']
    if debug:
        print(("--- (4) %s seconds ---" % (time.time() - start_time)))
    # Apply function to calculate mask value
#    df['mask'] = df.apply(mask_nighttime, axis=1)
    df['mask'] = df.apply(lambda x: mask_nighttime(
        x['lons'], x['lats']), axis=1)
    if debug:
        print(("--- (5) %s seconds ---" % (time.time() - start_time)))
    # Re-index by lat and lon
    df = pd.DataFrame(df['mask'].values, index=[df['lats'], df['lons']])
    if debug:
        print(("--- (6) %s seconds ---" % (time.time() - start_time)))
    # Unstack and return just as array
    df = df.unstack()
    marr = df.values
    if debug:
        print(("--- (end-7) %s seconds ---" % (time.time() - start_time)))

    return marr
