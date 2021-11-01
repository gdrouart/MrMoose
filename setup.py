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

reqs=['pathos>=0.2.8',
     'tqdm>=4.8.4',
     'scipy>=1.4.1',
     'numpy>=1.21.1',
     'emcee>=3.0.2',
     'corner>=2.2.1',
     'matplotlib>=3.4.2',
     'astropy>=4.2.1',
     'PyYAML>=5.4.1',
     'h5py>=3.3.0',
     'ultranest>=3.3.0',
     'cython>=0.29.24']

setup(name="MrMoose",
      version="2.0.0",
      author="Guillaume Drouart, Theresa Falkendal",
      author_emal="guillaume.drouart@curtin.edu.au, tfalkend@eso.org",
      description="MrMoose, multi-resolution, multi-source SED fitting code, with examples",
      license="GPUv3 Licence ",
      url="https://github.com/gdrouart/MrMoose",
      install_requires=reqs)
