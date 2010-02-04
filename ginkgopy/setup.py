#! /usr/bin/env python

###############################################################################
##
## GINKGO Biogeographical Evolution Simulator Post-Processing Library.
##
## Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
## This program is free software; you can redistribute it and#or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program. If not, see <http:##www.gnu.org#licenses#>.
##
###############################################################################

"""
Package setup and installation.
"""

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup
from setuptools import find_packages
from ginkgo import PACKAGE_VERSION

import sys
import os
import subprocess

script_names = []
setup(name='Ginkgo',
      version=PACKAGE_VERSION,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeet@ku.edu and mtholder@ku.edu',
      url='',
      description="""\
A library to faciliate setting up runs and processing results of the GINKGO Biogeographical Evolution Simulator""",
      license='GPL 3+',
      packages=['ginkgo'],
      package_dir={'ginkgo': 'ginkgo'},
      package_data={
        "" : ['doc/*'],
        "ginkgo" : ["tests/data/*"]
      },
      scripts = [('scripts/%s' % i) for i in script_names],
      test_suite = "ginkgo.tests",
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      long_description=open('README.txt', 'rU').read(),
      classifiers = [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics evolution biology biogeography',
      )
