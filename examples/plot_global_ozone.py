
from __future__ import absolute_import
import AC_tools as AC
import sys
sys.path=["/work/home/bn506/Python"]+sys.path
print sys.path

import pygchem




def main():

    # Specify the directory of the output data
    results_dir = "/work/home/bn506/temp/temp_shani_res"
    wd = results_dir

    # Get the GeosChem output
    #O3_data = AC.get_GC_output( results_dir, species='O3')
#    O3_data = AC.get_GC_output( results_dir, vars=['BIOGSRCE__ISOP'])
    O3_data = AC.get_GC_output( results_dir, vars=['vov'])


    # Get a 2d slice from the 3d array
    # Data is in the form lat,lon,alt
    O3_data = O3_data[:,:,0]

    print O3_data.shape
    O3_data = O3_data.T

    AC.map_plot( O3_data , res='0.5x0.625', wd=wd)
#    AC.map_plot( O3_data , wd=wd)

    import matplotlib.pyplot as plt
    plt.show()

if __name__ == '__main__':
    main()


