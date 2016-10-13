#!/usr/bin/env python

from distutils.core import setup
setup(  name='AC_tools',
        version='0.1',
        description='Atmospheric Chemistry tools',
        author='Tomas Shewen, Ben Newsome',
        author_email='Tomas.Sherwn@york.ac.uk',
        url='https://github.com/tsherwen/AC_tools',
        licence="why are they all so complicated - citations would be nice :-)",
        packages=['AC_tools'],
        install_requires=[
            'numpy',
            'wget>3.0',
            'filecmp',
            'pytest',
            'logging',
            'os',
        ],
        )

# Add any other install scripts here, like the user config setup that we want.


        
