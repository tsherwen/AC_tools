#!/usr/bin/env python
from setuptools import setup, find_packages
import os

VERSION = '0.1.1'
DISTNAME = 'AC_tools'
DESCRIPTION = "Atmospheric Chemistry (AC) tools"
AUTHOR = 'Tomas Sherwen, Ben Newsome'
AUTHOR_EMAIL = 'tomas.sherwen@york.ac.uk'
URL = 'https://github.com/tsherwen/AC_tools'
LICENSE = 'MIT'
PYTHON_REQUIRES = '>=3.5'

on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    INSTALL_REQUIRES = []
else:
    INSTALL_REQUIRES = [
                        'affine',
                        'cartopy',
                        'geopandas',
                        'matplotlib',
                        'netcdf4',
                        'numpy',
                        'pandas',
                        'pytest',
                        'rasterio',
                        'scipy',
                        'xarray'
                       ]

CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering',
]


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name=DISTNAME,
      version=VERSION,
      license=LICENSE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      classifiers=CLASSIFIERS,
      description=DESCRIPTION,
      long_description=readme(),
      python_requires=PYTHON_REQUIRES,
      install_requires=INSTALL_REQUIRES,
      url=URL,
      packages=find_packages())
