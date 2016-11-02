from ..funcs4GEOSC import *
import logging
import pytest

wd = '../data'

def test_get_surface_area():
    arr = get_surface_area(wd=wd)
    assert isinstance( arr, np.ndarray), 'Surface area not a numpy array.'
    assert (len(arr.shape)==3), 'Surface area does not have 2 dimensions.'


    # TO-DO
    # Test without passing a wd. Requires a data source. current implementation of this is a bit cludgey.
    arr = get_surface_area(updated=True)
    assert isinstance( arr, np.ndarray), 'Surface area not a numpy array.'
    assert (len(arr.shape)==3), 'Surface area does not have 2 dimensions.'
    
#    arr = get_surface_area(updated=True)
#    assert isinstance( arr, np.ndarray), 'Surface area not a numpy array.'
#    assert (len(arr.shape)==3), 'Surface area does not have 2 dimensions.'
    # Warning, this only works for tomas as requires setting up stuff.
    # The setup should not be needed if wd is passed about.

    return

def test_get_land_map():
    # Test specifying wd
    arr = get_land_map(wd=wd)
    assert isinstance( arr, np.ndarray), 'Land map not a numpy array.'
    assert (len(arr.shape) == 3), 'Land map has too many dimensions.'
    return

def test_get_O3_burden():
    var = get_O3_burden(wd=wd)
    assert (round(var.sum(),0) == round(376875.15625,0)), "The O3 burden is wrong ({var})".format(var=var)
    return

def test_get_air_mass_np():
    arr = get_air_mass_np(wd=wd)
    assert isinstance( arr, np.ndarray), 'Air mass array is not a numpy array'
    return

def test_get_GC_output():
    arr = get_GC_output(wd=wd, species='O3', category='IJ_AVG_S')
    assert isinstance( arr, np.ndarray), 'GC output is not a numpy array'
    assert round(arr.sum(),6)==round(0.14242639,6), "The ozone budget doesnt seem correct({bud})".format(bud=arr.sum())
    return

