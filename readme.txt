This ART data analyzer package is a python package containing a easy-to-use command line tool art_data that integrates with artn (http://normandmousseau.com/ART-nouveau-81.html) developed by Prof Mousseau and parallel computation, for the purpose of automating:

1) parallel computation of running Activation and Relaxation Techniques (ART) simulation jobs to generate ART data wrapped by python
2) post-processing ART data by filtering, calculating, visualizing, correlating various physical quantities change for large amount of event searches.

Currently, version 1.1 supports the automation of parallel ART simulations running and post-processing of ART data, archive and clean ART data for researchers to share data.
Now ART_data_analyzer support the visualization of potential energy landscape in terms of disconnectivity graph by integrating with pele package developed by David Wales group at Cambridge University (https://github.com/pele-python/pele). Future versions should implement the automation of integrating with LAMMPS, Ovito, pymatgen.

This python package is developed by Dr Liang Tian while he was originally at Dr Lin Li group at University of Alabama and now has been integrated into part of artn project in gitlab.
Developer: Dr Liang Tian.
Contact: liangtianisu@gmail.com


Cite this software:
Liang Tian, Lin Li, Jun Ding, Normand Mousseau, 
"ART_data_analyzer: automating parallelized computations to study the evolution of materials”, Software X, 9(2019), 238-243.

License:
This package is under GNU General Public License v3.0

Acknowledgement:

The author acknowledge many valuable discussions with Prof Mousseau and the financial support of U.S. DOE grant DOE: DE-SC0016164
The author acknowledge the feedback and suggestions from early stage users, e.g.Keshab Bashyal. 

Application and impact:
1) Calculate the activation energy barrier distribution for kinetics studies for a wide range of materials, 
e.g. crystalline materials, amorphous materials, crystalline defects (interfaces and grain boundaries)

2) correlating local structural characteristics with local properties from large amount of data to derive dynamical evolution of local structure-properties relation.


How to install artn and ART_data_analyzer?


A very detailed installation guide of artn is included in /ART_data_analyzer/ART_installation_guide/ along with environmental files (used to install and use artn) 
and various Makefiles for compiling lammps and artn.

Using the ART_data_analyzer package and its command line tool art_data need to correctly set up the python environment for your Operating system (though no need to compile and build python code).

This python package depends on some python packages such as numpy, pandas, matplotlib, python-tk,scipy, mpl_toolkits, scikit-learn, pathos, pyvoro, tess. 
First, user need to set up the python development environment for your OS, such as install python setuptools to get easy_install, pip. pip is the best way to handle additional python packages installation.

For Linux Ubuntu OS,
Install easy_install and pip by >>> sudo apt-get install python-pip python-dev build-essential
For Mac, install pip by >>> sudo easy_install pip

install previous mentioned python package by >>> python -m pip install --user numpy scipy matplotlib pandas scikit-learn pyvoro pathos tess. 

After installation of these dependencies, if there is no issue using art_data command line tool later, then it is good to go. 

However, sometimes it is common to see some errors due to the fact that some python packages, e.g. scikit-learn depends on another python package numpy and their version is not compatible. Solution: >>> pip install --upgrade --ignore-installed —-user numpy scipy matplotlib pandas scikit-learn pyvoro pathos tess
as seen in https://stackoverflow.com/questions/15418949/pip-wont-install-python-packages-locally-with-user . In this case, pip install --ignore-installed --user numpy will install a local user copy of numpy into $HOME/.local/lib/python2.7/site-packages, while ignoring the global copy of numpy located generally in /usr/lib/python2.7/site-packages. If for some specific cases, it needs a specific version of the package, we can do >>> pip install --user --force-reinstall numpy==1.15.4

If there is any error relating to any of the installed python packages mentioned above, check the bottom section on common issue while installing dependent python packages

If user want to use ART_data_analyzer to visualize the potential energy landscape by plotting disconnectivity graph, they should also install pele package by the following link: https://github.com/pele-python/pele. In Linux Ubuntu, following steps are needed:
1) pip install —-user networkx==2.2 sqlalchemy hungarian pyro4 brewer2mpl cython
2) git clone or download zip file of the repository from https://github.com/pele-python/pele
3) python setup.py build --fcompiler=gfortran
4) python setup.py install --user
After these steps, it should build fortran library for python package to use. The installed pele package will be located at $HOME/.local/lib/python2.7/site-packages/

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
	
	Based on my test, the sample size should be less than 225,000 atoms for pyvoro and tess to calculate voronoi indexes without being killed. However, I have successfully circumvent this intrinsic limitation of pyvoro and tess by avoiding calculating voronoi indexes for the whole sample. I take each of the interested points and find its neighbors in terms of a cube surrounding this point considering periodic boundary condition. Then I calculate the voronoi indexes of this small sample and only output the voronoi index of this point. This approach not only avoid the large sample memory error, but also improve the voronoi calculation speed for all states substantially. However, the regular nearest neighboring points surround the point in terms of sphere would not work due to the fact that spherical shape will cause the missing cell error as referred here:
https://github.com/joe-jordan/pyvoro/issues/18

Though sample size limitation has been addressed, it is still suggested that user should try their best to reduce the sample size since a large amount of events data with large sample size will occupy a large amount of disk space. For example, each configuration file of 310,000 atoms sample will take more than 20MB. A few thousand events will occupy 50-100 GB disk space.


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

Now support running ART simulations in slurm cluster environment by art_data -s input_sample_id.json --art --run —slurm
However, user need to coordinate with the slurm admistrator to set up the necessary environment for running artn before running in slurm environment.
In this case, the num_of_proc key in input_sample_id.json means the number of compute nodes
user requested to submit their sbatch jobs to. This will split all central atoms in central_atom_list key
into num_of_proc folds. Each fold will be allocated to one compute node as a single submitted sbatch job to use its all available cores to 
perform parallel computation using python multiprocesssing on this single compute node.

Some other errors may occur due to the usage of subprocess module to invoke sbatch job to run python multiprocessing module by multiprocessing.Pool class:
“””
"/share/apps/python/2.7/lib/python2.7/multiprocessing/pool.py", line 347, in _handle_results
    task = get()
TypeError: ('__init__() takes at least 3 arguments (1 given)', <class 'subprocess.CalledProcessError'>, ())
“””
This error is a python development bug according to this python issue post:
https://bugs.python.org/review/9400/diff/159/Lib/subprocess.py

User may request system administrator to modify the subprocess.py file /share/apps/python/2.7/lib/python2.7/subprocess.py by adding the following line
super(CalledProcessError, self).__init__(returncode, cmd, output)
 

If ART simulations are interrupted due to either human intervention or machine failure, user can check the current status of all tests
by —-check_tests_status, such as art_data -s input_sample_id.json --art --check_tests_status
This will create two input SETTINGS files called input_tests_done.json storing the id of finished tests and 
input_tests_undone.json storing the id of unfinished tests
	
For finished tests:
user can check the results of only finished tests as described in the
following post-processing section by e.g.
art_data -s input_tests_done.json --filter, art_data -s input_tests_done.json --eng --calc
	
For unfinished tests:
user can choose to delete these unfinished tests completely by art_data -s input_tests_undone.json --art --delete_tests
Then user can continue to run these unfinished tests from the beginning by art_data -s input_tests_undone.json --art --run

ART simulations are performed under constant volume so that the simulation box size does not change.

For ubuntu machine with small RAM (e.g.8GB), when running art in parallel in ubuntu box, some I/O error issues could occur during the writing of a single large configuration file (>20MB, >300,000 atoms) to obtain an incomplete/corrupted configuration file due to the insufficient RAM memory and swap space. 
It is suggested that the user should allow sufficient swap space according to their RAM as described in this article: https://www.cyberciti.biz/tips/linux-swap-space.html
For further information on swap (swap file or swap partition) in Ubuntu and how to increase the swap file size, check the following two articles:
https://help.ubuntu.com/community/SwapFaq
https://askubuntu.com/questions/927854/how-do-i-increase-the-size-of-swapfile-without-removing-it-in-the-terminal/927870

New features of artn:
1) local force evaluation:
LOCAL_FORCE .true.
will only perform the deformation and force evaluation inside an inner spherical region defined by INNER_REGION_RADIUS
with a shell of atoms surrounding this region to ensure forces acting on inner region are correct
Each configuration will be fully relaxed before writing into the configuration file to avoid building long range strain effect

2) remove the configuration files of rejected events to save disk space:
WRITE_REJECTED_EVENT:
default .true., all events configurations files are saved,
Max_num_events parameter now in bart.sh is total number of attempted events. 
 .false. it will overwrite the sad and min files for rejected events so that only saving the configuration files of accepted events.
Max_num_events parameter now in bart.sh is total number of accepted events. 
events.list file will always contain all accepted/rejected events so that you can see event history.


Post-processing tasks in a user workflow:

The following demonstrates the main workflow though more functionalities are available
1)filtering events:
	art_data -s input_sample_id.json --filter

2)perform activation and relaxation energy convergence statistical t tests:
	art_data -s input_sample_id.json --eng --calc
	
	art_data -s input_sample_id.json --eng --ttest_kfold OPTION k n
	
	art_data -s input_sample_id.json --eng --ttest OPTION PATH_1 PATH_2
	
	art_data -s input_sample_id.json --art --run_more N_TESTS
	
	art_data -s input_sample_id.json --update_input

After filtering and energy convergence tests, user can choose to delete unused ART data by
art_data -s input_sample_id.json —-art —-delete_unused_events
This will also update central_atom_list.json file automatically, save an original copy if necessary
art_data -s input_sample_id.json —-update_input
This will update the input_sample_id.json file automatically, save an original copy if necessary

3)run calculations and visualizations:
	 
	art_data -s input_sample_id.json --strain --calc
	art_data -s input_sample_id.json --strain -v
	art_data -s input_sample_id.json --strain --stats

4)determine the local atom indexes of each filtered event in list_of_test_id
by a machine learning outlier detection algorithm

5)find all local atom indexes for each of filtered events

6)run calculations for various user specified atom list 

For user’s convenience, it currently support local atoms by —-local, central atom by —-central, initial cluster atoms by —-initial,
max displaced atom by —-max_disp. 
Before using these arguments, user first need to find these various types of atoms by --find_local_index,--find_central_index,--find_triggered_cluster_atoms_index,--find_max_disp_index

For generality, it supports passing a list of atoms such as [1,6,8] as the value of atom_list key in input settings file


Speed demo:
ART running wrapper automatically implements the ART input files modification for each tests and parallelization of jobs submission. With the new feature of artn incorporating local force evaluation, a speed test of running 20,000 events and 100,000 events in a 10,000 atoms sample on a 24 cores 3.0 GHz, 64 GB RAM machine took about 12h and 3 days to accomplish, respectively. 
Usually, metallic glass sample needs about 20,000 initial events to be calculated to get sufficient data to reproduce the statistics.

Goal of future version:

The long term goal of this ART_data_analyzer package would be integrated into a bigger package that establish a true multiscale and multiphysics framework (with the option of being multi-stage machine learning to simplify complex physics details at each length scale for the purpose of minimizing computational cost for accelerating rational material design without losing the multi-scale physics picture) for at least the material design of metallic glass community, by integrating LAMMPS (already has a nice python wrapper), ART, ART_data_analyzer (python ART wrapper, post-processing with machine learning tool, maybe incorporating some features of data processing and visualization software Ovito by its python wrapper), Materials project python API pymatgen to query some materials properties data, 3D STZ dynamics. “Art is never finished, only abandoned.”


Scientific background:
The overall goal is to correlate the local change in atomic structures with the change of local properties (e.g. activation energy) by sampling sufficient number of local potential energy basins to reproduce the correct statistical distributions as training data for the machine learning model to train correlation models, by triggering various local events. By forcing to trigger a group of atoms surrounding central atom initially away from their minimum (1st step of ART perform minimization, after minimization, the configuration should be very close to the lammps input configuration to guarantee we are getting the activation energy barrier of initial lammps input sample), the ARTn 1st order saddle point convergence algorithm will climb the potential energy landscape to its nearest saddle state. During the initial to saddle process, the machine learning outlier detection algorithm will check how many atoms are inelastically rearranged during this trigger to estimate the number of involved atom(size of flow defect) during the activation of a potential hypothetical STZ (since we only call it STZ when it actually rearrange). Each triggering event corresponds to a potential individual STZ being activated to rearrange their atoms. Each activation of individual STZ is corresponding to a local beta relaxation. With these information probed by ART, thinking reversely, the physics is that the thermal energy kT relative to STZ activation energy distribution will determine the percentage of STZs will be activated. The percolation of activated STZs will transfer from local atomic rearrangement events, i.e. activation of individual STZs, into a global atomic rearrangement event associated with/characterized by a large global potential energy decrease, i.e. alpha relaxation, to finish the glass transition into their equilibrium liquid-like phase as we increase the thermal energy kT to activate more STZs. The critical temperature is the glass transition temperature. The deformation of glass is affected by this process since the external stress can not only directly change the rate of STZ activation in STZ dynamics, but also stress should be able to change the activation energy distribution probed by ART by changing the atomic structures, which is usually not incorporated. It is possible that stress may not change the average activation energy barrier though it can distort the potential energy landscape. It is also interesting to see softness (average rearrangement displacement magnitude of local involved atoms in STZs or nonaffine displacement D2min) correlated with local shear modulus, and activation energy.



This package has the following modules and sub-packages:

/src contains all the source code of various libraries, majorly in python, it may also contains other code such as Matlab or C code in the future

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
