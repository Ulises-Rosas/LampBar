import sys
import platform
import setuptools
from distutils.core import setup

if platform.architecture()[0] != '64bit':
    sys.stderr.write('Package requires a 64bit architecture')
    sys.stderr.flush()
    exit()

myos = sys.platform

if myos == 'darwin':
    bins = ['./ext_bin/darwin/agrep']

elif myos == 'linux' or myos == "linux2":
    bins = ['./ext_bin/linux/agrep']
    
else:
    sys.stderr.write('Package requires a UNIX-based OS')
    sys.stderr.flush()
    exit()

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="lampbar",
      version='0.1',
      long_description = readme,
      long_description_content_type='text/markdown',
      packages=['lampy'],
      data_files = [ ('bin', bins) ],
      classifiers=[
          'Programming Language :: Python :: 3'
      ]
    )


