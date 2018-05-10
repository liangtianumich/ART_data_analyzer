# set up environmental variables before using this package

# the directory path of the python package in user's machine
export MY_ART=~/Documents/alabama/ART_data_analyzer

# user need to specify the ART data directory
export DATA_DIR=/home/ltian/Documents/10^13/2000atoms/ART/Cu64Zr36_newpar/

# user need to specify the test directory for running tests
export TEST_DIR=$MY_ART/tests/

# add this package /src/python dir to PYTHONPATH so that the python modules can be searched
export PYTHONPATH=$MY_ART/src/python:$PYTHONPATH

# make the python running script as an exe file inside package dir, add this package dir to path
export PATH=$MY_ART/scripts:$PATH



