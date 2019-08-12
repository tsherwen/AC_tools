from ..plotting import *
import logging
import pytest
import os
logging.basicConfig(filename='test.log', level=logging.DEBUG)


wd = '../data'
out_dir = 'test_output'

slow = pytest.mark.skipif(
    not pytest.config.getoption("--slow"),
    reason="need --slow option to run"
)

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


@pytest.fixture()
def test_data():
    from ..GEOSChem_bpch import get_GC_output
    test_data = get_GC_output(wd, species='O3')
    return test_data


# Leave the default plot in as not slow.
def test_map_plot_default(test_data):
    print(test_data.shape)
    map_plot(test_data[:, :, 0, 0])
    return


@slow
def test_map_plot_transpose(test_data):
    map_plot(test_data[:, :, 0, 0].T)
    return


@slow
def test_map_plot_wd(test_data):
    map_plot(test_data[:, :, 0, 0], wd=wd)
    return


@slow
def test_map_plot_wrong_shape(test_data):
    with pytest.raises(AssertionError):
        map_plot(test_data[0, 0, :, :])
    return


@slow
def test_map_plot_none():
    with pytest.raises(AssertionError):
        map_plot(None)
    return


@slow
def test_save_plot_default(test_data):
    map_plot(test_data[:, :, 0, 0])
    save_plot()
    filename = "myplot.png"
    os.remove(filename)
    return


@slow
def test_save_plot_with_name(test_data):
    map_plot(test_data[:, :, 0, 0])
    save_plot(title="test_1")
    filename = "test_1.png"
    os.remove(filename)
    return


@slow
def test_save_plot_in_folder(test_data):
    map_plot(test_data[:, :, 0, 0])
    save_plot(title="test_3", location="new_folder")
    filename = os.path.join("new_folder", "test_3.png")
    os.remove(filename)
    os.rmdir("new_folder")
    return


@slow
def test_save_plot_with_filetypes(test_data):
    map_plot(test_data[:, :, 0, 0])
    save_plot(title="test_4", extensions=["pdf", "png"])
    filename = "test_4.png"
    filename2 = "test_4.pdf"
    os.remove(filename)
    os.remove(filename2)
    return

    # Test plot has been created, and then remove it"
#    filenames = [filename_0, filename_1, filename_2,
#                filename_3, filename_4, filename_5]
#    for filename in filenames:
#        assert os.path.isfile( filename ), "Failed to create {file}".format(file=filename)
#        os.remove(filename)
#
#    os.rmdir("new_folder")


pytest.main()
