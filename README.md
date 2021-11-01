
# MrMoose README - Current Version v2.0.0

## Table of contents:
* Presentation
* References
* Last updates
* Installation
* First run
* Known limitations
* Advanced features
* Compilation of tips and tricks


## Presentation:
MrMoose stands for Multi-Resolution Multi-Object/Origin Spectral Energy distribution fitting procedure.
In a nutshell, this code allows to fit user-defined models onto a set of multi-wavelength data using a
Bayesian framework. The code is shared with multiple examples to demonstrate its capabilities and allow
the user to adapt it more easily to their specific requirements. MrMoose is designed to permit user a
large freedom in a user-friendly fashion, but not being user-opaque. The code can therefore handle blended sources,
large variation in resolution, and even upper limits consistenly. It also generates a series of ouputs
allowing for an quick interpretation of the results. The code is using the emcee package, and saving the
emcee sampler object (.pkl) allowing users to tranfer the output to personnal graphical interface. Since v2.0, Ultranest, 
a much more efficient package for parameter space exploration is available. It is recommanded to prefer it over emcee
given the significant gain in performance and behaviour. 

## References: #
If using this code, please refer to the published paper or the ascl link of the code:
https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.4981D/abstract 
https://ui.adsabs.harvard.edu/abs/2018ascl.soft09015D/abstract

## Last updates, main modifications:
* v1.0.0: public release
* v1.0.1: correctif release
            - correction of the dump fit_struct in the .sav file to be human readable
            - addition of the model_struct containing parameters best fit values and uncertainties in .sav file
            - update of model library, correction of overflow problems in synchrotron laws
            - correction of graphic bug in displaying the walker chains plot (inverted color depending on AF_cut value)
* v1.1.0: update release
            - implementation of redshift as free parameters
            - correction of various bugs
            - homogeneisation of examples to free parameter redshift implementation
* v1.2.0: update to Python 3
            - change in file structure
            - add the conda environement file for easier installation
* v2.0.0: implementation of Ultranest sampler
            - refactoring of core functions for performance gain
            - v2.0 with Ultranest is about 100 times faster compared to v1.2 emcee

## Installation instruction:
Download MrMoose from github, either by running in a terminal (require git to be installed)
git clone https://github.com/gdrouart/MrMoose.git
(or as a standard .zip file from the same url in a browser)

Run the following in the directory, to install all required dependencies
python install setup.py

Add the path to MrMoose in your PYTHONPATH variable in your .bashrc or .profile to be called on demand
export PYTHONPATH=Your_path_to_MrMoose_directory:$PYTHONPATH

## First run
We summarise here the first run of MrMoose. Under ipython for instance
>run example1.py
will create three files, with the following path:
 - data/fake_source_ex1.dat
 - models/fake_source_ex1.mod
 - fake_source_ex1.fit

The .dat file contains all the data, the .mod file contains all the models (function calls) and
the interval to be considered for each parameter and the .fit file contains the setting to
perform the fit. 

To run your first fit, after launching ipython:
>import mrmoose
>tab = mrmoose.SED_fit('fake_source_ex1.fit')
will perform your first run and generate your first results! (several .pdf files and a .pkl file
in the outputs folder)

Make sure you also create a data, models and outputs folders if not already existing, otherwise MrMoose will crash!
>mkdir data
>mkdir models
>mkdir outputs

## Limitations, known bugs:
- Initial values set as the median of the interval of parameters - not true anymore with Ultranest
- Require to underestimate the parameters(especially normalisation factor) in case of the presence of upper limits
- Parallelisation on one source only is not effective, but working efficiently for sample (one source per core)
- Only works in a python 2.7 environment and emcee 2.1, it is recommanded to create a virtual 
    environment to avoid conflicts with the most recent package versions - solved in v1.2.0

## Future developments:
- Implementation of checkpoints to re-start chain convergence in case of crash/stops
- Allowing different prior on the parameter (only uninformative, uniform prior in v1.0)
- Implement Jupyter interface
- Transform MrMoose in a package
- Implement use of the logging system
- Implement template libraries of non-linear models to be fit along
- Migration to Python 3 - Version 1.2
- Allowing different redshift for different components, and allowing redshift as a free parameter - added in Version 1.1
- Move the advanced feature to the setting file (.fit) as optional parameters

## Advanced features
- AF_cut: in the mrmoose.py file, the option "AF_cut" can be modified to filter manually walkers below a
given acceptance fraction(AF) value. The default value is -1, where the code filter automatically chains
with AF<mean(AF)/2, but is not necessarily optimal in certain cases. To obtain a view of the distribution
of AF values, change the "histo" value to True (located just below the AF_cut). 
- histo: This will enable the plot of a histrogram of the AF values during the convergence plot procedure. 
- layout: in the mrmoose.py the layout option allow for different style: presentation and talk. Customisation
is possible by adding your own keyword and options - following the rcParams dictionary of matplotlib format -
in the mm_utilities.py file in the function named "format_plot_output".
- AICc(in development): a number to compare different model combination directly. See Akaike Information Criteria
definition for a description of this diagnostic tool for model comparison. 

## Tips and tricks:
We provide here some extra tips to manipulate ouputs:

To compress the pdf into more handy version:
 - pdf to compressed pdf
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH -sOutputFile=output.pdf input.pdf
 - pdf to png (or jpeg), just play around with the density and resizing to obtain a good quality image much lighter! 
convert -density 200 -resize 200% input.pdf -quality 100 output.png
