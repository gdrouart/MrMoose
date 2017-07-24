###############################################
# ########################################### #
# # MrMoose README - Current Version v1.0.0 # #
# ########################################### #
###############################################

Table of contents:
* Presentation
* References
* Last updates
* Installation
* First run
* Known limitations
* Advanced features
* Compilation of tips and tricks


#################
# Presentation: #
#################
MrMoose stands for Multi-Resolution Multi-Object/Origin Spectral Energy distribution fitting procedure.
In a nutshell, this code allows to fit user-defined models onto a set of multi-wavelength data using a
Bayesian framework. The code is shared with multiple examples to demonstrate its capabilities and allow
the user to adapt it more easily to their specific requirements. MrMoose is designed to permit user a
large freedom in a user-friendly fashion, but not being user-opaque. The code can therefore handle blended sources,
large variation in resolution, and even upper limits consistenly. It also generates a series of ouputs
allowing for an quick interpretation of the results. The code is using the emcee package, and saving the
emcee sampler object (.pkl) allowing users to tranfer the output to personnal graphical interface. 


###############
# References: #
###############
If using this code, please refer to the published paper of the code:
arxiv pdf: 
bibtex: 


#####################################
# Last updates, main modifications: #
#####################################
* v1.0.0: public release


#############################
# Installation instruction: #
#############################
Download the MrMoose package and run the following in the directory:
python install setup.py

This should install all required dependencies. After this, MrMoose can be imported as
a package and used in this directory, using classical python commands


#############
# First run #
#############
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


############################
# Limitations, known bugs: #
############################
- Initial values set as the median of the interval of parameters
- Require to underestimate the parameters(especially normalisation factor) in case of the presence of upper limits
- Parallelisation on one source only is not effective, but working efficiently for sample (one source per core)


########################
# Future developments: #
########################
- Implementation of checkpoints to re-start chain convergence in case of crash/stops
- Allowing different prior on the parameter (only uninformative, uniform prior in v1.0.0)
- Implement template libraries of non-linear models to be fit along
- Migration to Python 3
- Allowing different redshift for different components, and allowing redshift as a free parameter
- Move the advanced feature to the setting file (.fit) as optional parameters

#####################
# Advanced features #
#####################
- AF_cut: in the mrmoose.py file, the option "AF_cut" can be modified to filter manually walkers below a
given acceptance fraction(AF) value. The default value is -1, where the code filter automatically chains
with AF<mean(AF)/2, but is not necessarily optimal in certain cases. To obtain a view of the distribution
of AF values, change the "histo" value to True (located just below the AF_cut). 
- histo: This will enable the plot of a histrogram of the AF values during the convergence plot procedure. 
- layout: in the mrmoose.py the layout option allow for different style: presentation and talk. Customisation
is possible by adding your own keyword and options - following the rcParams dictionary of matplotlib format -
in the mm_utilities.py file in the function named "format_plot_output".

####################
# Tips and tricks: #
####################
We provide here some extra tips to manipulate ouputs:

To compress the pdf into more handy version:
 - pdf to compressed pdf
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH -sOutputFile=output.pdf input.pdf
 - pdf to png (or jpeg)
convert -density 200 -resize 200% input.pdf -quality 100 output.png
