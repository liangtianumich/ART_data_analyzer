This ART data analyzer package is a python package that automates 1) parallel computation of running ART simulation jobs to generate ART data as a python wrapper 2) data analysis of various physical quantities change between two atomic configurations of each possible event search of using Activation and Relaxation Techniques (ART). The installation guide of ART is included in /ART_installation_guide along with environmental files when installing ART and various Makefiles for compiling lammps and ART. Currently, version 1.0 supports the ART running wrapper and post-processing of ART calculated data.  The automation implementation of integrating with LAMMPS, Ovito, pymatgen will come in future versions.

ART running wrapper automatically implements the ART input files modification for each tests and parallelization of jobs submission. With the new feature incorporating local force evaluation in new ART,a speed test of running 100,000 events on a 24 cores 3.0 GHz, 64 GB RAM machine took only about 3 days to accomplish. Usually, metallic glass sample needs about 20,000 initial events to be calculated to get sufficient data to reproduce the statistics. After the user finished the following steps to correctly set up the environment for this package, user can use the art_data —-desc to read the details about how to use this package, for example, check the convergence of data and run more data if necessary. 

The overall goal is to correlate the local change in atomic structures with the change of activation energy by sampling sufficient number of local potential energy basins to reproduce the correct statistical distributions as training data for the machine learning model to train correlation models, by triggering various local events. By forcing to trigger a atom initially away from their minimum (1st step of ART perform minimization, after minimization, the configuration should be very close to the lammps input configuration to guarantee we are getting the activation energy barrier of initial lammps input sample), the ART algorithm will climb the potential energy landscape to its nearest saddle state and check how many atoms are rearranged during this trigger to estimate the number of involved atoms/size of STZ during the activation of a potential hypothetical STZ (since we only call it STZ when it actually rearrange, a potential STZ with too high activation energy should not be called STZ). Each triggering event corresponds to a potential individual STZ being activated to rearrange their atoms. Each activation of individual STZ is corresponding to a local beta relaxation. With these information probed by ART, thinking reversely, the physics is that the thermal energy kT relative to STZ activation energy distribution will determine the percentage of STZs will be activated. The percolation of activated STZs will transfer from local atomic rearrangement events, i.e. activation of individual STZs, into a global atomic rearrangement event associated with/characterized by a large global potential energy decrease, alpha relaxation, finish the glass transition into their equilibrium liquid-like phase as we increase the thermal energy kT to activate more STZs. The critical temperature is the glass transition temperature. The deformation of glass is affected by this process since the external stress can not only directly change the rate of STZ activation in STZ dynamics, but also stress should be able to change the activation energy distribution probed by ART by changing the atomic structures, which is usually not incorporated. It is possible that stress may not change the average activation energy barrier though it can distort the potential energy landscape. It is also interesting to see softness (average rearrangement displacement magnitude of local involved atoms in STZs or nonaffine displacement D2min) correlated with local shear modulus, and activation energy 


Currently, it can calculate and visualize the information such as atomic displacement, atomic strain (e.g. the volumetric and shear strain), voronoi index as a statistical distribution for all filtered events. It will also implement the physically meaningful atomic structure descriptor for representing the atomic structure within local environment to be correlated with the local properties.

This package has the following modules and sub-packages:

/src contains all the source code of various libraries, majorly in python, it may also contains other code such as Matlab or C code

data_reader: extracting data from various data file formats (default format lammps dump file), read data into a pandas.Dataframe

event_selector: load the accepted events and select events passed stage 1 two criteria first to save these events into selected_events.json 
Selected events will be further filtered by the redundant event search function, which implement 
redundancy check of all possible pairs of events that has already passed stage 1 two criterias

Util: utilities containing various classes and functions for other modules

calculator package:
	
	strain_calculator: calculate the atomic strains for all events or a user customized subset of events

	voronoi analysis: calculate the voronoi index and voronoi class for all atoms for all events for a suer customized subset of events

	stress_calculator,  coming soon
		implemented by LAMMPS based on “How thermally activated deformation starts in metallic glass”
	
	shear_modulus calculator, coming soon 
		implemented by LAMMPS based on “Configurational dependence of elastic modulus of metallic glass”
	
	D2min calculator:
		non-affine displacement modified in terms of “dynamic” in state rather than dynamic in time
	
	/correlation_model:
		contains modules to perform correlation analysis between these local atomic physical quantities derived from local atomic structures, 
		and local atomic properties such as activation energy of local potential energy barrier basin
	
	/local_atom_env_calculator: calculate the local environment quantities for each atom, such as the fundamental quantities bond distance, bond orientation, coordination number to replace a complete set of structural functions represented by the radial basis functions and spherical harmonics, that will reduce the risk of overfitting due to the limited number of training data, simplify the physics by machine learning models for computational costs during later correlation analysis
	
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


/examples: contains examples of executable running scripts and demo data. The demo example system is CuZr metallic glass.

/scripts: contains executable scripts that will be directly invoked on the command line

/tests: contains unittests to test various modules



Operating system:
This python package depends on some python packages such as numpy, pandas, matplotlib, python-tk,scipy, mpl_toolkits, scikit-learn, pathos, pyvoro. User need to set up the python development environment for your OS, such as install python setuptools to get easy_install, pip. Pip is the best way to handle additional python packages installation.

For Mac or Linux Ubuntu OS,
Install easy_install and pip by sudo apt-get install python-pip python-dev build-essential, install previous mentioned python package by python -m pip install --user numpy scipy matplotlib jupyter pandas scikit-learn pyvoro pathos. Mac O.S mpl_toolkits under matplotlib library may not be a package and need user to manually make it a package to import. User need to verify the successful installation of these dependencies in their OS. The package has been tested on python 2.7. Python 3 user may need further tests, especially on print(), xrange to range, os.makedirs(path_to_dir,exist_ok=True) to handle race condition when in event level parallelization: https://stackoverflow.com/questions/12468022/python-fileexists-error-when-making-directory or https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist. Pyvoro released on PyPI is version 1.3.2. https://pypi.org/project/pyvoro/, which does not support python3. The user need to install from source by python setup.py install --user on this branch https://github.com/joe-jordan/pyvoro/tree/feature/python3 , which support python3 https://github.com/joe-jordan/pyvoro/pull/11. For python3 user, they may need to refer to another python interface to voro++ library Tess to support python3: https://github.com/wackywendell/tess

A few notes:
if pip install does not work, especially for pip version 9.0.1, user need to upgrade their pip pip-10.0.1 according to here: https://stackoverflow.com/questions/49768770/not-able-to-install-python-packages-ssl-tlsv1-alert-protocol-version.  Then it can install package such as scikit-learn by pip install —-user scikit-learn


User also need to install pathos package for a unified event level parallel computation by using the util.operation_on_events function by pip install --user pathos. Pathos can serialize any function that are not defined at the top level of the module that greatly favor the implementation. Check: https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error. It seems that any file that call util.operations_on_events() will need to import os, pickle numpy as np.

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


For user who actively use other python and packages version and need to create an isolated python package environment just for using this ART_data analyzer package without interfering with your global python environment, please use virtualenv path/to/ART_data_analyzer to create a new virtual python environmental for the path/to/ART_data_analyzer (if need inheriting python from global environment python, use --system-site-packages flag like virtualenv --system-site-packages path/to/ART_data_analyzer), then it will install /bin /include /lib under path/to/ART_data_analyzer to include python packages for this virtual environment. The bin/activate is bash script used to activate this virtual python environment by source /ART_data_analyzer/bin/activate.  This will change the $PATH (add /path/to/ART_data_analyzer/bin) and $VIRTUAL_ENV = /path/to/ART_data_analyzer to let you use the virtual python environment. And type deactivate to end this virtual environment. Create a alias start_ART = “source /ART_data_analyzer/bin/activate” in .bash_profile or .profile to activate faster. More details check: https://virtualenv.pypa.io/en/stable/userguide/
https://packaging.python.org/tutorials/installing-packages/


How to install/use this python package:
This package is written in python scripting language to analyze the atomic configuration data generated by ART atomistic simulation software. There is no need to compile and build python code. This package can be put in any directory. Ensure to source the environment.sh before using this package. The purpose is to create the necessary environmental variables PYTHONPATH/PATH/TEST_PROJ for current bash sessions to find the python packages/exe scripts/test directory. The user need to change these environmental variables to point to their correct ART_data_analyzer package and the dir to ART data on their machine

After the user source environment.sh, the user can use the art_data exe to parse the command line arguments of user input to run all functions of this package. 
To see a list of options, type in terminal art_data -h,
It will prompt all available optional arguments that can use for various functions of the python package to analyze the art_data.
The user can use art_data —-desc to see how to run various optional arguments.

More options will be implemented in art_data, in order to achieve the following goal.

“Art is never finished, only abandoned.”

The long term goal of this ART_data_analyzer package would be integrated into a bigger package that establish a true multiscale and multiphysics framework (with the option of being multi-stage machine learning to simplify complex physics details at each length scale for the purpose of minimizing computational cost for accelerating rational material design without losing the multi-scale physics picture) for at least the material design of metallic glass community, by integrating LAMMPS (already has a nice python wrapper), ART (which need a python wrapper), ART_data_analyzer(data processing by machine learning, maybe incorporating some features of data processing and visualization software Ovito by its python wrapper), Materials project python API pymatgen to query some materials properties data, 3D STZ dynamics.


The following sections is written in early development stage to record some nice design principles and can be skipped if necessary

How to use the package as a user work flow?
Though we say this is a data analysis package, the data here simply means ART raw data. To  obtain useful information from these atomistic simulation raw data, we need to perform the following steps:

After you set up the environment variables to get correct path to data and ART_data_analyzer.

First (data preprocessing), run event_filter.py to select non-redudant events saved in final_selected_events.json, or run your calculation, e.g.strain_calc.py if the user think that event redundancy check in event_filter.py would not filter many events.

Second(data calculation), if you run your calculation in first step, now you need to run event_filter to get correct statistical distribution. mp module will parallelize the calculation in a single test one by one

Third(data visualization and data analysis), run your data visualization (strain visualization) or data analysis module (correlation model) to extract the correct information, all these modules have been implemented with pathos.multiprocessing module through the util.operation_on_events function, where operation is a function to do data visualization or data analysis on a single event, pathos.mp will parallelize this operation function on single event one by one. This can be different from test level parallelization if test contains more than one final selected events. Correlation model module can be used to determine the average number of local involved atoms from initial to saddle configuration along with their histograms by implementing sklearn.linear_model.RANSACRegressor or sklearn.svm.LinearSVR to find outliers, for the user to get reproducible result from RANSAC, user need to set up the seed integer in random_state option of RANSAC class and most importantly, choose the residual_threshold: maximum residual for a data sample to be classified as an inlier. By default the threshold is chosen as the MAD (median absolute deviation) of the target values y. Currently, the user is suggested to use LinearSVR to find the average number of local involved atoms since in this fitted linear model of LinearSVR, average linear slope is constant that best reflect the physical constraint. Another two functions are available to get the fixed number of outliers using sklearn.ensemble.IsolationForest and sklearn.neighbors.NearestNeighbors.

Fourth, after the user get the average number of local atoms, then the user can involve the local atoms strain calculation mode for all events by letting atom_list to be a dictionary with key “local”. This will perform the strain calculation for all local atoms defined as the number of atoms that are closest to the triggered atom which are defined in ART input file bash.sh. The strain_calculation.py module has implemented the finding of these local atoms automatically. If all atomic strains are saved into previous strain_results.pkl and displacement_results.pkl, it will load the strains of local atoms and save them into local_strain_results.pkl and local_displacement_results.pkl. If no any calculations exists, it will calculate and save into the local_strain_results.pkl and local_displacement_results.pkl.

Then the user can repeat the data visualization and data analysis (to develop correlation model) for locally involved atoms to avoid the local effect being averaged out by a large amount of atoms. To obtain the accurate value of average number of locally involved atoms, the user should consider use the relative residual_threshold parameter sweep for locally involved atoms for all final selected events to find a consistent stopping criteria to decide the ave number of local atoms for this specific material system, as a standard work flow. Currently, the stopping criteria is two consecutive relative residual threshold data points whose slope is less than -15, i.e.,0.15 local atom number change/0.01 relative residual-(y_max-y_min). It is possible that this stopping criteria (e.g. will also be used to determine the ave number of local atoms for data calculated under different initial condition, such as the different initial sample generated by LAMMPS MD under different cooling rate.


Executable files:


The following exe scripts are the individual scripts to run individual functions at the early development stage. The user only need to use art_data to run various functions.

Currently, strain_calc.py is an exe python file, it performs displacement and atomic strain calculations and save all results and plots automatically in their corresponding locations. The calculation results have been rigorously verified with the results by Ovito.

voronoi_index_calc.py calculates the voronoi index for user specified events. The geometric voronoi diagram calculation under periodic boundary condition is implemented by the pyvoro package: https://github.com/joe-jordan/pyvoro. It is available in PyPI and can be easily installed by pip and easy to use. 

event_filter.py ran rigorous redundancy check for each possible pair events which has passed the two stage 1 criteria, this script can ran before strain calculation to calculate minimal number of meaningful events or after strain calculation to help select the filtered events for the correct statistical distribution, 
whether to use this script depends on the percentage of events filtered by this criteria. Generally, the user is encouraged to use this script on a subset of tests such as 10 tests to check how many events are filtered. If negligible, It is possible to ignore this criteria.

strain_visualization.py generate the plots of all available calculated events that passed two stage 1 criteria, even when there is unfinished calculations or calculations are going on since it ignored uncalculated events, user can customize a subset of their interested tests to extract statistics

Before running this script, after user specify all the environmental variables in environmental.sh to match their own machine and source it. 

Nice features: 
support parallel calculation by multiprocessing module, the user can specify the number of processes in the input file

all calculations in the fly of strain calculations will be saved into a pkl file, such as nn_results.pkl, strain_results.pkl, and displacement_results.pkl for initial to saddle configuration, and saddle to final configuration for each event. When rerunning this calculation, these files will be directly read to prevent redundant running. It will skip the tests whose log.file or configuration files (e.g either min1000.dump or min1000) does not exist. It will also skip the test who do not have either accepted or selected events after calculation.

A re_calc argument is included with default False to use the existing calculations. However, when user want to try calculate with a new set of input parameters, they may need to invoke the re_calc to True.

All plots for visualizing the statistics of various physical quantities are automatically plotted and saved into their corresponding event locations that are specified by event_initial_str_sad_str_fin_str: for example, event_min1001_sad1002_min1002

These plots are automatically updated without redoing the calculation so that user can modify the plotting function to customize their own plot styles

The mean, std, max of displacement, shear strain, volumetric strain of each event are plotted into a histogram for all events, these plots are also automatically updated.

All events are accepted/selected based on their individual modules and have minimal interaction with the calculation modules. User can customize their own event selection criteria.

The user can specify the cut-off distance cut-off, the simulation box dimension box_dim and the total number of tests, num_of_procs in their data directory
{'cut_off':cut_off_distance,'box_dim':box_dim,'num_of_tests':num_of_tests,’num_of_procs’:num_of_procs}

If the num_of_procs == 1, it will run calculation in single process mode

If num_of_procs >1, it will invoke the multiprocessing module Pool class
num_of_procs can be found by nproc --all or grep -c processor /proc/cpuinfo in linux

In Mac, sysctl -n hw.ncpu to get logical CPUs

Crossplatform option here, default is to use all CPU cores by multiprocess.cpu_count() as implemented as default option here.

An example of cut_off_distance dict: {(1,1):3.7,(1,2):3.7,(2,2):3.7}
Means atom_type_1 and atom_type 1 cut_off_distance is 3.7, type 1 and type 2 is also 3.7 etc

Box_dim only support orthogonal currently, more to be implemented in the future

User need to input the number of tests in their data directory


This input for running the data calculation will be saved into os.environ['TEST_DIR']
as a pickle file for future reference if needed


 



