# set up environmental variables before using this package

# the directory path of the python package in user's machine
export MY_ART=$HOME/Documents/alabama/ART_data_analyzer

# user need to specify the ART data directory
export DATA_DIR=$MY_ART/examples/

# user need to specify the test directory for running tests
export TEST_DIR=$MY_ART/tests/

# No need to change, add this package /src/python dir to PYTHONPATH so that the python modules can be searched
export PYTHONPATH=$MY_ART/src/python:$PYTHONPATH

# No need to change, make the python running script as an exe file inside package dir, add this package dir to path
export PATH=$MY_ART/scripts:$PATH


# for running ART
# set up the path to ART executable
export ART_EXE=$HOME/artn/source/ARTn_exec

export ART_INPUT=$DATA_DIR



