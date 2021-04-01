'''
NUPACK setup
'''

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from setuptools.dist import Distribution

from setuptools.command.install import install
# To use a consistent encoding
from codecs import open

# https://stackoverflow.com/questions/60379212/how-to-create-a-full-wheel-with-abi-tag
class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name"""
    def has_ext_modules(foo):
        return True

class InstallPlatlib(install):
    def finalize_options(self):
        install.finalize_options(self)
        if self.distribution.has_ext_modules():
            self.install_lib = self.install_platlib

setup(
    name='nupack',
    version='@PROJECT_VERSION@',
    description='Nucleic Acid Package',
    url='www.nupack.org',
    package_data={'nupack': ['cpp.so', 'parameters/*']},
    packages=find_packages(include=('nupack**',)),
    scripts=[],
    distclass=BinaryDistribution,
    cmdclass={'install': InstallPlatlib},
    install_requires=[
        'scipy>=1.0',
        'numpy>=1.17',
        'pandas>=1.1.0',
        'jinja2>=2.0',
    ]
)

