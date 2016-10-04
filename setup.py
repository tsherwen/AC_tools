"""
Setup script to install everything needed.
Downloads the data needed for the tests.
"""

def get_test_files():                                                                       
    """                                                                                    
    Downloads all the test dataset files using rsync.                                      
    """                                                                                    
    import os                                                                              
    if not os.path.isfile("test_files/ctm.nc"):                                            
        os.system("""(cd tests/test_files && wget -r -nH -nd -np -R --no-parent --reject "index.html*" http://atmosviz1.york.ac.uk/~bn506/data/GC_funcs_test/test_files/)""")
                                                           

if __name__ == "__main__":
    get_test_files()

    
