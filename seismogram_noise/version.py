#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 2
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev0'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 4 - Beta",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "seismogram_noise: Create synthetic noise time series"
# Long description will go up on the pypi page
long_description = """
Seismogram noise
========

License
=======
``SEISMOGRAM_NOISE`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.
All trademarks referenced herein are property of their respective holders.
Copyright (c) 2018 - Simon Staehler (mail@simonstaehler.com)
"""

NAME = "seismogram_noise"
MAINTAINER = "Simon Staehler"
MAINTAINER_EMAIL = "mail@simonstaehler.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/sstaehler/seismogram_noise"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Simon Staehler"
AUTHOR_EMAIL = "mail@simonstaehler.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['seismogram_noise']
PACKAGE_DATA = {'seismogram_noise': [pjoin('data', '*')]}
REQUIRES = ['numpy', 'scipy']
