
"""
Test suit for plotting.py
"""

import logging
from .. import plotting

def test_plot_settings():
    plot_settings = plotting.plot_settings()
    assert type(plot_settings.tick_size)  == int, "Tick size is not an integer"
    assert type(plot_settings.title_size) == int, "Title size is not an integer"
    assert type(plot_settings.filename)   == str, "Output_file is not a string"

    return

def test_multi_line_plot():
    
    data_1 = {
        'name'      : 'Line 1',
        'x_data'    : [1,2,3,4,5],
        'y_data'    : [2,2,2,2,2],
        }
    data_2 = {
        'name'      : 'Line 2',
        'x_data'    : [1,2,4,5],
        'y_data'    : [0,1,3,4],
        }
    data = [data_1, data_2]
    x_ticks = [[1,2,3,4,5],["a","b","c","d","e"]]
    plot_settings = plotting.plot_settings()
    plot_settings.filename = "test_outputs/multi_line_plot.png"
    plotting.multi_line_plot( data, x_ticks=x_ticks, plot_settings=plot_settings)

    return
