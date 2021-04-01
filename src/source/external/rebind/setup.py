import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="rebind",
    version="0.0.0",
    author="Mark Fornace",
    description=("Python bound C++ unit testing"),
    license="Not decided",
    keywords="C++ unit test framework",
    url="http://packages.python.org/an_example_pypi_project",
    packages=['rebind'],
    long_description=read('README.md'),
)
