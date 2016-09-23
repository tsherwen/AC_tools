import pytest
def pytest_addoption(parser):
    parser.addoption("--remake_ctm", action="store_true",
        help="remake ctm.nc tests")
