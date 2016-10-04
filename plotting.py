"""
Settings and functions to help with plotting.
"""
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from beeprint import pp


class plot_settings:
    '''
    Class containing plot settings which can be passed to plotting functions
    default populated values:
    tick_size, title_size, filename, tick_axis, tick_type, labelleft
    
    '''
    def __init__(self):
        '''
        Creates the default plot settings to use in matplotlib
        '''

        logging.info('Creating a plot settings class')

        self.tick_size  = 18
        self.title_size = 26
        self.filename   = "test.png"
        self.tick_axis  = "y"
        self.tick_type  = "major"
        self.labelleft  = "no"
        self.zeroline   = None
        self.legend_location = None
        self.title = None
        self.x_label = None
        self.y_label = None
        self.x_ticks = None
        self.y_ticks = None
        
        return


def setup_plot( plt,  plot_settings ):

    if plot_settings==None:
        from .plotting import plot_settings as get_plot_settings
        plot_settings = get_plot_settings()
    # rename for faster coding
    p_sets = plot_settings


    
    if p_sets.zeroline==True:
        plt.axhline(y=0, color='k')

    if not p_sets.x_ticks==None:
        plt.xticks( x_ticks[0], x_ticks[1] )
    if not p_sets.y_ticks==None:
        plt.yticks( y_ticks[0], y_ticks[1] )
    if p_sets.legend_location==None:
        pass
    elif p_sets.legend_location=='right':
        plt.legend(loc="center left", bbox_to_anchor=(1,0.5))
    else:
        plt.legend(loc=p_sets.legend_location)

    if not p_sets.title==None:
        plt.title( p_sets.title, size=p_sets.title_size )

    if not p_sets.x_label == None:
        plt.xlabel( p_sets.x_label )
    if not p_sets.y_label == None:
        plt.ylabel( p_sets.y_label )


    return 

def filled_multi_line_plot( dataset, x_ticks=None, y_ticks=None, plot_settings=None 
                    ):
    '''
    Creates a set of line plots from a list of lists.
    '''

    logging.info('Plotting a multi line plot.')

    # Check if we have plottable data.
    assert type(dataset) == list, "The x_data is not a list!"
    for data in dataset:
        assert type(data) == dict, "An item in x_data is not a dict."

    if plot_settings == None:
        from .plotting import plot_settings as get_plot_settings
        plot_settings = get_plot_settings()

    plot_settings.x_ticks = x_ticks
    plot_settings.y_ticks = y_ticks

    figure = plt.figure()


    ls = '-'
    label = None

    first_run = True
    rainbow = ['r','c','m','y','k','b','g']*10
    for data, color in zip(dataset,rainbow):
        x_data = data['x_data']
        y_data = data['y_data']
        if 'name' in data:
            label = data['name']
        if 'ls' in data:
            ls = data['ls']


        if first_run:
            first_run=False
            old_data = data
            plt.plot( x_data, y_data, label=label, color=color, ls=ls )
            continue

        y_data_upper = data['y_data']
        y_data_lower = old_data['y_data']

        plt.plot( x_data, y_data, label=label, color=color )
        plt.fill_between( x_data, y_data_lower, y_data_upper,
                         facecolor=color, label=label )
        old_data = data



    
    configure_plot( plt, plot_settings )

    figure.savefig(plot_settings.filename, bbox_inches='tight')
    print "figure saved to {filename}".format(filename=plot_settings.filename)
    logging.info("figure saved to {filename}".format(filename=plot_settings.filename))

    return

def num2col( number ):
    '''
    Get a color from a number
    '''

    colors = ['r','g','b','m','c','y','k', 'brown', 'orange']
    number = number % len(colors) 

    color = colors[number]
    return color


def multi_line_plot( dataset, x_ticks=None, y_ticks=None, plot_settings=None 
                    ):
    '''
    Creates a set of line plots from a list of lists.
    '''

    logging.info('Plotting a multi line plot.')

    # Check if we have plottable data.
    assert type(dataset) == list, "The x_data is not a list!"
    for data in dataset:
        assert type(data) == dict, "An item in x_data is not a dict."


    # Set up the plot.
    figure = plt.figure()


    # Plot the data
    label=None
    ls='-'
    lc=None

    for data in dataset:
        x_data = data['x_data']
        y_data = data['y_data']
        if 'name' in data:
            label = data['name']
        if 'ls' in data:
            ls = data['ls']
        if 'lc' in data:
            lc = data['lc']
            if isinstance(lc, int):
                lc = num2col(lc)
        plt.plot( x_data, y_data, label=label, ls=ls, color=lc )


    # Save the figure

    setup_plot( plt, plot_settings )
    figure.savefig(plot_settings.filename, bbox_inches='tight')
    print "figure saved to {filename}".format(filename=plot_settings.filename)
    logging.info("figure saved to {filename}".format(filename=plot_settings.filename))

    return

