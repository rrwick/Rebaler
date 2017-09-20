# -*- coding: utf-8 -*-

"""
Run 'python3 setup.py install' to install Rebaler.
"""

import sys
import re
from setuptools import setup


# Make sure this is being run with Python 3.4 or later.
if sys.version_info.major != 3 or sys.version_info.minor < 4:
    sys.exit('Error: you must execute setup.py using Python 3.4 or later')

version = re.search('^__version__\s*=\s*"(.*)"', open('rebaler/rebaler.py').read(), re.M).group(1)
short_descr = re.search("ArgumentParser\(description='(.*)'", open('rebaler/rebaler.py').read(), re.M).group(1)

with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')

setup(name='Rebaler',
      description=short_descr,
      long_description=long_descr,
      url='http://github.com/rrwick/rebaler',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      version=version,
      license='GPL',
      packages=['rebaler'],
      install_requires=['biopython'],
      entry_points={'console_scripts': ['rebaler = rebaler.rebaler:main']})
