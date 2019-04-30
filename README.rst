AC_tools - Atmospheric Chemistry (AC) tools
======================================

This module contains a functions and scripts used for 
working with atmospheric model output and observational data. 
Many functions are included for working with global and regional 
chemical transport model (CTM) ouput from the GEOS-Chem model.

This package started as a just collection of scripts that were
found to be useful for work in atmospheric chemistry and now
simply aims to contain functionality outside the remit of the 
more specialised community packages (e.g PyGChem_, xbpch_, and 
gcpy_) and use the existing Python stack (e.g. dask_, xarray_, 
pandas_). `Pull Requests are 
welcome! <https://github.com/tsherwen/AC_tools/pulls>`_

Installation
------------

**AC_Tools** is currently only installable from source. To do this, you
can either clone the source directory and manually install::

    $ git clone https://github.com/tsherwen/AC_tools.git
    $ cd AC_tools
    $ python setup.py install

or, you can install via pip directly from git::

    $ pip install git+https://github.com/tsherwen/AC_tools.git

Quick Start
-----------

Functions within **AC_Tools** can be used for various tasks for handling model output and observations. 

An exmample would be importing NetCDF files or converting ctm.bpch files from a directory of GEOS-Chem_ output (with ``tracerinfo.dat`` and ``diaginfo.dat`` files). 


.. code:: python

    import AC_tools as AC
    folder = '<folder containing GEOS-Chem output>'
    # Get the atmospheric ozone burden in Gg O3 as a np.array
    array = AC.get_O3_burden(folder)
    print( "The ozone burden is: {burden}".format(burden=array.sum()))
    # Get surface area for resolution 
    s_area = get_surface_area(res)[..., 0]  # m2 land map
    # Get global average surface CO 
    spec = 'CO'
    array = AC.get_GC_output(wd=folder, vars=['IJ_AVG_S__{}'.format(spec)])
    ratio = AC.get_2D_arr_weighted_by_X(array, res='4x5', s_area=s_area) 
    print( "The global average surface mixing ratio of {spec} (ppbv) is: {ratio}".format(spec=spec, ratio=ratio*1E9))
    
    
Usage
------------

Example analysis code for using AC_tools is available in the 
scripts folder. 

For more infomation, please visit the AC_tools_wiki_.


License
-------

Copyright (c) 2015 `Tomas Sherwen`_

This work is licensed under a permissive MIT License.

Contact
-------

`Tomas Sherwen`_ - tomas.sherwen@york.ac.uk

.. _`Tomas Sherwen`: http://github.com/tsherwen
.. _conda: http://conda.pydata.org/docs/
.. _dask: http://dask.pydata.org/
.. _licensed: LICENSE
.. _GEOS-Chem: http://www.geos-chem.org
.. _xarray: http://xarray.pydata.org/
.. _pandas: https://pandas.pydata.org/
.. _gcpy: https://github.com/geoschem/gcpy
.. _PyGChem: https://github.com/benbovy/PyGChem
.. _xbpch: https://github.com/darothen/xbpch
.. _AC_tools_wiki: https://github.com/tsherwen/AC_tools/wiki
