"""
MrMoose: A multi-wavelength, multi-resolution, multi-sources fitting code.
Allow to simultaneously fit multiple models in order to interpret the spectral
energy distribution of an astrophysical source from simple to more complex case, 
making use of a bayesian framework.

Copyright (C) 2017 Guillaume Drouart, Theresa Falkendal

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import os
import sys
from setuptools import setup

# to update
reqs=['pathos>=0.2.0',
     'tqdm>=4.8.4',
     'scipy>=0.14.0',
     'guppy>=0.1.10',
     'numpy>=1.9.1',
     'emcee>=2.2.1',
     'corner>=2.0.1',
     'pycallgraph>=1.0.1',
     'matplotlib>=1.5.1',
     'astropy>=2.0rc1',
     'PyYAML>=3.12']

setup(name="MrMoose",
      version="1.2.0b",
      author="Guillaume Drouart, Theresa Falkendal",
      author_emal="guillaume.drouart@curtin.edu.au, tfalkend@eso.org",
      description="MrMoose, multi-resolution, multi-source SED fitting code, with examples",
      license="GPUv3 Licence ",
      url="https://github.com/gdrouart/MrMoose",
      install_requires=reqs)
