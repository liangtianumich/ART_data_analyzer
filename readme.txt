This ART data analyzer package is a python package containing a easy-to-use command line tool art_data that integrates with artn (http://normandmousseau.com/ART-nouveau-81.html) developed by Prof Mousseau and parallel computation, for the purpose of automating:

1) parallel computation of running Activation and Relaxation Techniques (ART) simulation jobs to generate ART data wrapped by python
2) post-processing ART data by filtering, calculating, visualizing, correlating various physical quantities change for large amount of event searches.

Currently, version 0.0.1 supports the automation of parallel ART simulations running and post-processing of ART data. Future versions should implement the automation of integrating with LAMMPS, Ovito, pymatgen.

This python package is originally developed by Dr Lin Li group at University of Alabama and now has been integrated into part of artn project in gitlab.
Lead developer: Dr Liang Tian.
Contact: liangtianisu@gmail.com

Acknowledgement:

The author acknowledge many valuable discussions with Prof Mousseau and the financial support of U.S. DOE grant DOE: DE-SC0016164

Application and impact:
1) Calculate the activation energy barrier distribution for kinetics studies for a wide range of materials, 
e.g. crystalline materials, amorphous materials, crystalline defects (interfaces and grain boundaries)

2) correlating local structural characteristics with local properties from large amount of data to derive dynamical evolution of local structure-properties relation.


How to install artn and ART_data_analyzer?


A very detailed installation guide of artn is included in /ART_data_analyzer/ART_installation_guide/ along with environmental files (used to install and use artn) 
and various Makefiles for compiling lammps and artn.

Using the ART_data_analyzer package and its command line tool art_data need to correctly set up the python environment for your Operating system (though no need to compile and build python code).

This python package depends on some python packages such as numpy, pandas, matplotlib, python-tk,scipy, mpl_toolkits, scikit-learn, pathos, pyvoro. 
First, user need to set up the python development environment for your OS, such as install python setuptools to get easy_install, pip. pip is the best way to handle additional python packages installation.

For Mac or Linux Ubuntu OS,
Install easy_install and pip by >>> sudo apt-get install python-pip python-dev build-essential

install previous mentioned python package by >>> python -m pip install --user numpy scipy matplotlib jupyter pandas scikit-learn pyvoro pathos. 

After installation of these dependencies, if there is no issue using art_data command line tool later, then it is good to go.

If there is any error relating to any of the installed python packages mentioned above, check the bottom section on common issue while installing dependent python packages

For Windows O.S user, check the bottom section on setting up python environment for Windows users.



How to use ART_data_analyzer python package/art_data command line tool?

Before using this package, make sure that you either source /ART_data_analyzer/ART_installation_guide/art_install.sh or copy the content in art_install.sh into .profile (need to log-out and re-log-in for ubuntu since ubuntu default bash is not a log-in shell). This is because the content in art_install.sh is not only used during installation of artn, but also needed when using artn and lammps to find the openmpi and openblas library path.

This package can be put in any directory. Currently, for the convenience of integrating with /artn, /ART_data_analyzer is located under $HOME/artn

1) modify /ART_data_analyzer/environment.sh:
	modify DATA_DIR to your desired directory path to save the ART data.
	
	modify MY_ART if necessary, ensure that /ART_data_analyzer is in $MY_ART.
	e.g. export MY_ART=$HOME/artn/ART_data_analyzer

	Generally, we do not need to modify other environmental variables,
	can be the same as default. Current default is that ART_INPUT storing path_to_input_files 
	are set up to be the same as path_to_data_dir so that all ART input files
	should be located under path_to_data_dir.
	
	The default sample is lammps_data type of sample, whose name must be
	modified to conf.lammps. This lammps data conf.lammps file can be
	converted from dump file by Ovito easily.
	export ART_SAMPLE=conf.lammps
	export SAMPLE_TYPE=lammps_data
	
	Sometimes, if user need to use the dump file sample, they can modify
	the ART_SAMPLE to the dump file name (e.g.conf.dump),  SAMPLE_TYPE=dump
	export ART_SAMPLE=conf.dump
	export SAMPLE_TYPE=dump

2) source the /ART_data_analyzer/environment.sh before using this package. The purpose is to create the necessary environmental variables PYTHONPATH/PATH/ for current bash sessions to find the python packages/exe scripts.

3) Now user can use the command line tool art_data to parse the command line arguments of user input to run all functions of this package. 

To see a list of options, type in terminal >>> art_data -h
It will prompt all available optional arguments that can be used for various functions of the python package.

The user can type >>> art_data —-desc 
to see a detailed guide on how to use various optional arguments.

Below are some demos of art_data can do, always check a complete guide by art_data —-desc:

input_sample_id.json is an example input SETTINGS file

automate running ART in parallel:
	art_data -s input_sample_id.json --art --run

Post-processing tasks in a user workflow:
1)filtering events:
	art_data -s input_sample_id.json --filter
2)perform activation and relaxation energy convergence statistical t tests:
	art_data -s input_sample_id.json --eng --calc
	
	art_data -s input_sample_id.json --eng --ttest_kfold OPTION k n
	
	art_data -s input_sample_id.json --eng --ttest OPTION PATH_1 PATH_2
	
	art_data -s input_sample_id.json --art --run_more N_TESTS
	
	art_data -s input_sample_id.json --update_input

3)run calculations and visualizations:
	 
	art_data -s input_sample_id.json --strain --calc
	art_data -s input_sample_id.json --strain -v
	art_data -s input_sample_id.json --strain --stats

4)determine the local atom indexes of each filtered event in list_of_test_id
by a machine learning outlier detection algorithm

5)find all local atom indexes for each of filtered events

6)run calculations for local atoms only


Speed demo:
ART running wrapper automatically implements the ART input files modification for each tests and parallelization of jobs submission. With the new feature of artn incorporating local force evaluation, a speed test of running 20,000 events and 100,000 events in a 10,000 atoms sample on a 24 cores 3.0 GHz, 64 GB RAM machine took about 12h and 3 days to accomplish, respectively. 
Usually, metallic glass sample needs about 20,000 initial events to be calculated to get sufficient data to reproduce the statistics.

Goal of future version:

The long term goal of this ART_data_analyzer package would be integrated into a bigger package that establish a true multiscale and multiphysics framework (with the option of being multi-stage machine learning to simplify complex physics details at each length scale for the purpose of minimizing computational cost for accelerating rational material design without losing the multi-scale physics picture) for at least the material design of metallic glass community, by integrating LAMMPS (already has a nice python wrapper), ART, ART_data_analyzer (python ART wrapper, post-processing with machine learning tool, maybe incorporating some features of data processing and visualization software Ovito by its python wrapper), Materials project python API pymatgen to query some materials properties data, 3D STZ dynamics. “Art is never finished, only abandoned.”


Scientific objective:
The overall goal is to correlate the local change in atomic structures with the change of activation energy by sampling sufficient number of local potential energy basins to reproduce the correct statistical distributions as training data for the machine learning model to train correlation models, by triggering various local events. By forcing to trigger a group of atoms surrounding central atom initially away from their minimum (1st step of ART perform minimization, after minimization, the configuration should be very close to the lammps input configuration to guarantee we are getting the activation energy barrier of initial lammps input sample), the ART algorithm will climb the potential energy landscape to its nearest saddle state and check how many atoms are rearranged during this trigger to estimate the number of involved atoms/size of STZ during the activation of a potential hypothetical STZ (since we only call it STZ when it actually rearrange). Each triggering event corresponds to a potential individual STZ being activated to rearrange their atoms. Each activation of individual STZ is corresponding to a local beta relaxation. With these information probed by ART, thinking reversely, the physics is that the thermal energy kT relative to STZ activation energy distribution will determine the percentage of STZs will be activated. The percolation of activated STZs will transfer from local atomic rearrangement events, i.e. activation of individual STZs, into a global atomic rearrangement event associated with/characterized by a large global potential energy decrease, alpha relaxation, finish the glass transition into their equilibrium liquid-like phase as we increase the thermal energy kT to activate more STZs. The critical temperature is the glass transition temperature. The deformation of glass is affected by this process since the external stress can not only directly change the rate of STZ activation in STZ dynamics, but also stress should be able to change the activation energy distribution probed by ART by changing the atomic structures, which is usually not incorporated. It is possible that stress may not change the average activation energy barrier though it can distort the potential energy landscape. It is also interesting to see softness (average rearrangement displacement magnitude of local involved atoms in STZs or nonaffine displacement D2min) correlated with local shear modulus, and activation energy.



This package has the following modules and sub-packages:

/src contains all the source code of various libraries, majorly in python, it may also contains other code such as Matlab or C code

data_reader: extracting data from various data file formats

event_selector: filtering events module, implementing stage 1 two criteria and stage 2 one criteria


Util: utilities containing various classes and functions for other modules

calculator package:
	event_energy_calculator: calculate the event activation and relaxation energy, perform energy convergence statistical t test. 
	
	strain_calculator: calculate the atomic strains and displacement

	voronoi analysis: calculate, classify and plot the voronoi index

	stress_calculator,  coming soon
		implemented by LAMMPS based on “How thermally activated deformation starts in metallic glass”
	
	shear_modulus calculator, coming soon 
		implemented by LAMMPS based on “Configurational dependence of elastic modulus of metallic glass”
	
	D2min calculator:
		non-affine displacement modified in terms of “dynamic” in state rather than dynamic in time

ART_wrapper package:
	run: automate art input files modification and parallel computation of ART jobs
	
correlation_model package:
	contains modules to perform correlation analysis between these local atomic physical quantities derived from local atomic structures,
	and local atomic properties such as activation energy of local potential energy barrier basin
	
	/local_atom_env_calculator: calculate the local environment quantities for each atom, such as the fundamental quantities bond distance, bond orientation, coordination number to replace a complete 	set of structural functions represented by the radial basis functions and spherical harmonics, that will reduce the risk of overfitting due to the limited number of training data, simplify the 	physics by machine learning models for computational costs during later correlation analysis
	
	bond_distance_calculator: calculate all bond distance for each atom
		an extra function will read the results file and extract statistics for correlation analysis
	bond_orientation_calculator: calculate the bond orientation for all three atom pair for each center local atom
		an extra function will read the results file and extract statistics
	coordination_number_calculator: calculate and summarize the number of different atomic species around the center atom
		an extra function will read the results file and extract statistics	
	
	calculate local atomic physical quantities that are derived from local atomic structure/environment
		
visualizer package:
	general_visualizer: functions for other visualizer modules to use 

	strain_visualizer: visualize the atomic strain related quantities
	
	voronoi_visualizer: visualize the voronoi cell classification


/examples: contains examples of demo data. The demo example system is CuZr metallic glass.

/scripts: contains executable scripts that will be directly invoked on the command line

/tests: contains unittests to test various modules



Common issues while installing dependent python packages:

if pip install does not work, especially for pip version 9.0.1, user need to upgrade their pip pip-10.0.1 according to here: https://stackoverflow.com/questions/49768770/not-able-to-install-python-packages-ssl-tlsv1-alert-protocol-version.  Then it can install package such as scikit-learn by pip install —-user scikit-learn

Mac O.S mpl_toolkits under matplotlib library may not be a package and need user to manually make it a package to import. User need to verify the successful installation of these dependencies in their OS. 

The package has been tested on python 2.7. Python 3 user may need further tests, especially on print(), xrange to range, os.makedirs(path_to_dir,exist_ok=True) to handle race condition when in event level parallelization: https://stackoverflow.com/questions/12468022/python-fileexists-error-when-making-directory or https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist. 

Pyvoro released on PyPI is version 1.3.2. https://pypi.org/project/pyvoro/, which does not support python3. The user need to install from source by python setup.py install --user on this branch https://github.com/joe-jordan/pyvoro/tree/feature/python3 , which support python3 https://github.com/joe-jordan/pyvoro/pull/11. For python3 user, they may need to refer to another python interface to voro++ library Tess to support python3: https://github.com/wackywendell/tess

User also need to install pathos package for a unified event level parallel computation by using the util.operation_on_events function by pip install --user pathos. Pathos can serialize any function that are not defined at the top level of the module that greatly favor the implementation. Check: https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error. It seems that any file that call util.operations_on_events() will need to import os, pickle numpy as np.

It seems that when installing the latest pathos, it will install a new version of dill, which will cause some import error. This can be fixed by using dill version 0.2.7.1 by pip uninstall dill and then pip install —-user dill==0.2.7.1


Setting up python environment for Windows users:

For windows 10 user, you can search Ubuntu app in Microsoft App Store to use Ubuntu app like regular linux Ubuntu. Make sure you are running Windows Powershell as admistrator by right click and copy commands https://docs.microsoft.com/en-us/windows/wsl/install-win10
If you encounter error like wslregisterdisribution failed with error ox800700b7
First, try search windows feature, turn off windows subsystem for linux (WSL), reboot, then turn on windows subsystem for linux, reboot. This should work. If not work, try reinstall ubuntu app from app store, see https://github.com/Microsoft/WSL/issues/2982
Now your home directory will be /home/unix_user_name, where your windows C drive can be accessed by /mnt/c/.

For windows 7, Cygwin should be installed to use bash.

For windows Cygwin user, set up python development environment to run python script in cygwin can be done by installing necessary python packages by downloading the official Cygwin setup_x86_64.exe (64 bit)or setup_x86.exe (32 bit) from the official Cygwin website and choosing the related packages. setup.exe performs like an official “package manager” for Cygwin, through there is no true package manager in Cygwin. For example, running the setup.exe to install the necessary packages needed for python development environment are 
\All\Devel\make
\All\Devel\gcc-core
\All\Devel\gcc-fortran
\All\Devel\gcc-g++
\All\Python\python2
\All\Python\python2-pip
Additional python packages may be needed to install some packages, such as installing pandas library may need the following necessary packages:
https://stackoverflow.com/questions/34341112/pandas-build-on-cygwin
python2-devel
numpy
python2-six
python2-wheel
python2-setuptools
python2-pip
python2-cython
wget
There is a package called apt-cyg acting like the package manager apt-get in Ubuntu that might be useful: https://superuser.com/questions/304541/how-to-install-new-packages-on-cygwin


Isolated python environments for different version of python:

For user who actively use other python and packages version and need to create an isolated python package environment just for using this ART_data analyzer package without interfering with your global python environment, please use virtualenv path/to/ART_data_analyzer to create a new virtual python environmental for the path/to/ART_data_analyzer (if need inheriting python from global environment python, use --system-site-packages flag like virtualenv --system-site-packages path/to/ART_data_analyzer), then it will install /bin /include /lib under path/to/ART_data_analyzer to include python packages for this virtual environment. The bin/activate is bash script used to activate this virtual python environment by source /ART_data_analyzer/bin/activate.  This will change the $PATH (add /path/to/ART_data_analyzer/bin) and $VIRTUAL_ENV = /path/to/ART_data_analyzer to let you use the virtual python environment. And type deactivate to end this virtual environment. Create a alias start_ART = “source /ART_data_analyzer/bin/activate” in .bash_profile or .profile to activate faster. More details check: https://virtualenv.pypa.io/en/stable/userguide/
https://packaging.python.org/tutorials/installing-packages/
