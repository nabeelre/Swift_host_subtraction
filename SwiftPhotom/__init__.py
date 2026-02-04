#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Giacomo Terreran (2021)
#
# This file is part of Swift_host_subtraction
#
# Swift_host_subtraction is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Swift_host_subtraction is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Swift_host_subtraction.  If not, see <http://www.gnu.org/licenses/>

"""Swift_host_subtraction
"""

from ._version import get_versions

__version__ = get_versions()['version']
__author__ = 'Giacomo Terreran <gterreran@lco.global>'
__credits__ = ['Peter Brown <grbpeter@yahoo.com>']

# Re-export commonly used submodules and functions
from . import commands, help, uvot
from .download import get_swift_data, download_swift_data, create_run_files
from .photom_host import run_swift_photom

__all__ = [
    'commands',
    'help',
    'uvot',
    'get_swift_data',
    'download_swift_data',
    'create_run_files',
    'run_swift_photom',
]

del get_versions
