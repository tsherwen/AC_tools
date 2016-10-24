import AC_tools as AC

# Specify the working directory
wd = "example_data"

# Get the GeosChem species data from the wd
my_data = AC.get_GC_output( wd, species='O3')

# Get a 2d slice from the 3d array
O3_data = O3_data[:,:,0]

# Create the plot
AC.map_plot( O3_data)

# Save the plot and show it.
AC.save_plot("my_plot")
AC.show_plot()


