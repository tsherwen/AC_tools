#####
# This file contains the settings used for py.test
#####

import pytest
import logging
import os
import urllib2

#FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"            
FORMAT = "%(filename)s:%(lineno)s - %(funcName)s() : %(message)s"
test_file_dir = 'test_files'                                                    

logging.basicConfig(filename='test.log',level=logging.DEBUG, format=FORMAT)
logging.getLogger().setLevel(logging.DEBUG)                                     

                                             


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true",
        help="remake ctm.nc tests")

def pytest_configure():

    # Make sure we are in the correct folder for the test.
    dirname = os.path.split(os.getcwd())[1]
    if not dirname=='Tests':
        pytest.exit("Not running in the Tests folder!")
                                                  
    """                                                                                    
    Downloads all the test dataset files using rsync.                                      
    """                                                                         
                                                                                
    test_files = ['test.nc', 'test.bpch','tracerinfo.dat','diaginfo.dat']       
                                                                                
    url_base =  'http://atmosviz1.york.ac.uk/~bn506/data/AC_tools/'             
    test_file_dir = 'test_files'                                                
                                                                                
                                                                                
    if not os.path.exists(test_file_dir):                                       
        os.makedirs(test_file_dir)                                              
                                                                                
    for file_name in test_files:                                                
        file_path = os.path.join(test_file_dir, file_name)                      
        if not os.path.isfile(file_path):                                       
            my_file = open(file_path, 'wb')                                     
            logging.debug(file_name + " not found. Downloading now.")           
            url = url_base + file_name                                          
            file_data = urllib2.urlopen( url ).read()                           
            my_file.write(file_data)                                            
            my_file.close()                                                     
                                                                                
            logging.debug(file_name + " downloaded.")                           
                                                                                
    return           
