# -*- test-case-name: p3d -*-
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License
# See LICENSE for more details.

"""
Python Protein Framework.
"""

__all__ = ["protein","vector","atom","tree","library","geo"]

# Ensure the user is running the version of python we require.
import sys

if not hasattr(sys, "version_info") or sys.version_info < (2,6):
    raise RuntimeError("p3d requires Python 3.0 or later.")
del sys

# Imports
from p3d import *
