import sys
from os.path import dirname, abspath, join

# Make sure that the application source directory (this directory's parent) is
# on sys.path.

here = join(dirname(dirname(abspath(__file__))), 'lib')
sys.path.insert(0, here)
