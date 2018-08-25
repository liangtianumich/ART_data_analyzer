1. Download
ART source code can be downloaded from Prof Normand’s website:
http://normandmousseau.com/ART-nouveau-81.html
Choose the most recent stable version “ART nouveau linked to the LAMMPS library”

For developers who want to contribute to ART, the development of ART is available in gitlab private repository so that developer can contact Prof Normand at nm@normandmousseau.com for access to the newest version of ART. 

After you get the source code tar.gz or zip file, extract it, move it to any location you want, here we suggested move it to $HOME/ART_stable,  where $HOME is /home/$USER in ubuntu and Mac. 

2 Installation of all dependent libraries: 

The installation of ART follows the section III DOC/guide.pdf (version March 19 2018 or early), however, there are some things need special attention and need corrections in this guide:

On Linux Ubuntu:

We will compile ART with gnu compiler: gcc, g++, gfortran, make sure you have all these compiler installed and included in your default environmental variable PATH by typing in the command line prompt: which gcc, which g++, which gfortran. If the command line returned a path to these executable, you are good to proceed. Otherwise, you have to install gcc, g++ and gfortran. 

The author can also choose to use the intel compiler icc, icpc, ifort as described in Doc/guide.pdf, but the author need to install the Intel compiler set by themselves since it is not installed by default. The vendor compiler like intel compiler may perform better, and is generally available in cluster by load intel compiler modules. 

The ART installation depends on the following libraries, whose installation sequences are important to a successful installation of ART. In clusters, these software libraries are available by load corresponding modules so that no need to install by yourself.

With some modifications from IIIB of guide.pdf, we are building the ART-depended libraries on Ubuntu as follows:

First, openmpi :
Your computer may already have an openmpi, to make sure you are using the openmpi we are going to install,  it is recommended to install this openmpi to a unique path, such as $HOME/installed_openmpi

Download the most recent stable openmpi from https://www.open-mpi.org

Extract it. 

Pre-build it by ./configure -prefix=$HOME/installed_openmpi CC=gcc CXX=g++ F77=gfortran FC=gfortran 
Here, it is important to add the CC CXX F77 FC option to make sure openmpi is wrapping around the right compiler (avoid the error of No underlying complier was specified in the wrapper compiler data file, mpicc-wrapper-data.txt, which is located in $HOME/installed_openmpi/share/openmpi) and also ensure that mpi.mod file will be present after final installation (it could raise an error can not use mpi can not find mpi.mod error when not using these option).

After this command finishes (if everything went well), install openmpi by typing make all install

Installation will be complete after about 20 minutes with a lot of screen output.

You can experiment with the make program parallel building -j option to invoke multiple core building option, however, make -j option will fail to build ART since the order of the dependencies to build ART can be important as will be described later.


Second, openblas, 
Download openblas source code from https://www.openblas.net/

Extract it.

cd into the main dir of extracted openblas (make sure gfortran is installed, or it will not compile lapack library and give error undefined reference to these functions dgeev and dsysv etc), type
make FC=gfortran 
After a lot of screen output, follow the bottom instruction (make PREFIX=/path/to/your/installation install), here we will install openblas to a unique path so that it will not interfere with your system openblas library that could be located in /usr/lib/

make PREFIX=$HOME/installed_openblas install
It will install openblas in $HOME/installed_openblas

Third(optional), fftw, the user do not have to install fftw3 since lammps has a built in KISS FFT library (see comments of lammps Makefile with fftw option), but if user want to install fftw by themselves (since fftw3 is faster than publically-available FFT implementations), they can do in the following:

Download source code from http://www.fftw.org/download.html
Extract it
Cd into main dir of extracted fftw,
./configure –prefix=$HOME/installed_fftw
make
make install
It will install fftw

3. Set up PATH and LD_LIBRARY_PATH environmental variables
art_install.sh is attached along with this dir to set up the PATH and LD_LIBRARY_PATH
In ubuntu, please pay attention LD_LIBRARY_PATH, not LD_LYBRARY_PATH in guide.pdf
In Mac, change LD_LIBRARY_PATH to be DYLD_LIBRARY_PATH, not DYLD_LYBRARY_PATH in guide.pdf

If you want to use your own installed fftw, uncomment lines in art_install.sh

In the terminal, type source art_install.sh

This will set up the PATH and LD_LIBRARY_PATH for the shell to find the exe or library we just installed in step 2. 

Make sure we are doing step 3 before compiling lammps library to ensure we are using the same openmpi for both exe and library header in both lammps and ART installation.
This can be double checked by which mpif90, or mpicxx, mpic++, mpicc, mpirun, which should give the same path as the installed path $HOME/installed_openmpi/bin

4 building lammps library
Since we will use jpeg and png dump image function in lammps, we will build lammps with
Jpeg and png library. First we need to install these two libraries by:

sudo apt-get install libjpeg-dev libpng-dev
If it says it can not find the package name, then do this first: sudo apt-get update

Download lammps source code
Extract it. 
Cd into /src, type make command to check a complete list of options for make command.
make is an building automation exe program default installed in /usr/bin, which takes the main lammps_dir/src/Makefile to check the options of make command in current lammps/src dir

User can check which packages will be included in the compiled lammps library by make package-ps
Then, if user want to install extra packages, first check if this package need external library by typing make package. It will show whether this package is standard or user package, whether this package need system libraries, or lammps included libraries or external libraries that user has to install for themselves. It would take extra efforts to install packages with external libraries, such as VORONOI. In the /src dir, type make lib-voronoi to see the instructions on how to download, build, install external libraries or go to the lib/voronoi/README to check the instructions. You can check the instruction of how to install external libraries for packages on lammps website: http://lammps.sandia.gov/doc/Section_start.html#start-3-3

Here for the basic needs and a quick installation, it is suggested to load the standard package without any library by make yes-standard; make no-lib
Some packages are the basic need in order for the lammps to run. 

If the user forget to add the desired package for lammps to run lammps input file so that ART can not read some commands of in.lammps, they need to re-compile lammps with those additional packages and also has to re-compile ART.

After you load all needed packages and their dependent libraries, you can start building lammps library by:

Copy the attached Makefile.lammps_for_art_ubuntu into the lammps_dir/src/MAKE/MINE

Note: 
1) if you want to use your installed fftw, uncomment line about 
FFT_INC =    
FFT_PATH = 
FFT_LIB =

2) JPG_INC is path to jpeg.h or png.h, whose default installation location is in /usr/include in ubuntu, change if necessary  
JPG_PATH is path to libjpeg.a or libpng.a, whose default location is /usr/lib/x86_64-linux-gnu in ubuntu, change if necessary

Cd into lammps_dir/src dir, type make mode=lib lammps_for_art_ubuntu

This will build lammps library with all needed packaged and a liblammps_lammps_for_art_ubuntu.a file should exist in lammps_dir/src folder


5. Build ART with lammps library:

Copy the Makefile.ubuntu_lammps in this dir into the ART_dir/source/MAKE
Cd into ART_dir/source, type make ubuntu_lammps
This should work. If not, the LAMMPS_USER_PATH is not set up correctly. Try either modify Makefile.ubuntu_lammps or Copy the liblammps_lammps_for_art_ubuntu.a and liblammps.a that links to liblammps_lammps_for_art_ubuntu.a into the ART_dir/source
Note:
1)Change the lammps_dir name if necessary, here lammps dir lammps-16Mar18 installed in $HOME
 LAMMPS_USER_PATH = $HOME/lammps-16Mar18/src

2) if you want to use your installed fftw, uncomment line about 
FFT_INC =    
FFT_PATH = 
FFT_LIB =

3)
JPG_INC is path to jpeg.h or png.h, whose default installation location is in /usr/include in ubuntu, change if necessary
  
JPG_PATH is path to libjpeg.a or libpng.a, whose default location is /usr/lib/x86_64-linux-gnu in ubuntu, change if necessary

Make -j option should not be used here, since it will compile the ART dependencies in the wrong order. Also no need to use -j since it compiles very quickly 

The art should be built successfully with exe ARTn_exec that links to exe ARTn_ubuntu_lammps. The final line of terminal output should look like
====
ln -sf ./ARTn_ubuntu_lammps_v1672 ./ARTn_exec
====

If there is error, please double check if the steps are followed successfully, especially about the installation of openmpi, openblas etc.


6. Run ART:

For ART to run, we need to have tcsh rather than sh, check if it exists by which tcsh, if it returns a path like /usr/bin/tcsh or /bin/tcsh, then both csh and tcsh is installed. 

If it does not exists, install both csh and tcsh by sudo apt-get install tcsh


Note for running ART:
The user need to source art_install.sh every time to set up the library path to openmpi and openblas whenever they want to USE the ART and/or lammps. 

It is suggested that the user can copy the content in art_install.sh into .profile in linux ubuntu or .bash_profile in Mac to make the environment set-up for ART and lammps permanent, may need log out current user account and log in to make it effective for ubuntu since the default ubuntu bash is not a log-in shell by default while Mac bash is a log-in shell by default. Whether log-in or non log-in of current shell (saved in $SHELL) can be checked by echo $0. Log-in shell will return -bash, non-log in return bash. If $SHELL is not bash, it is suggested to change the default shell to bash by chsh -s /bin/bash . It is not suggested to put these non-bash specific commands into bash specific .bashrc which is executed for any non-login shell. Other types of shell besides bash such as sh, ksh can also read non-bash specific commands in .profile. .profile usually also source .bashrc. If .bash_profile exists, .profile will not be read.



 







