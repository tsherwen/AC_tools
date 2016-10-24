import AC_tools as AC

# Download the example data if it is not already downloaded.
from AC_tools.Scripts import get_example_files

# Add some logging
import logging
#logging.basicConfig(filename='example.log',level=0)     


# Specify the working directory
wd = "example_data"

# Get the GeosChem species data from the wd
my_data = AC.get_GC_output( wd, species='O3')

# Get a 2d slice from the 3d array
my_data = my_data[:,:,0]

# Create the plot
AC.map_plot( my_data)

# Save the plot and show it.
AC.save_plot("my_plot")
AC.show_plot()


