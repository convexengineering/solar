"""Standard Python setup script for Solar Repo"""
from __future__ import print_function
import sys
from distutils.core import setup

LONG_DESCRIPTION = """
Solar is a repo designed for sizing analysis for solar, long-endurance aircraft.
"""

LICENSE = """MIT License

Copyright (c) 2017 Convex Engineering

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. """

setup(
    name="solar",
    description="Geometric programming model of solar long endurance aircraft",
    author="MIT Convex Engineering",
    author_email="mjburton@mit.edu",
    url="https://www.github.com/convexengineering/solar",
    install_requires=["numpy", "scipy", "gpkit", "pandas", "gpfit", "gpkitmodels"],
    version="0.0.0.0",
    packages=["solar"],
    package_data={"solar": ["*.csv"]},
    license=LICENSE,
    long_description=LONG_DESCRIPTION,
)
