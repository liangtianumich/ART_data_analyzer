#!/usr/bin/env python
import os
import pickle
import json
import pandas as pd
import matplotlib.pyplot as plt
from calculator.strain_calculator import strain_calculator_run_all_tests

path_to_data_dir = os.environ['DATA_DIR']

# save some backup files that might be useful after running the tests
path_to_test_dir = os.environ['TEST_DIR']

if not os.path.exists(path_to_test_dir):
	os.makedirs(path_to_test_dir)

# make this an input file to load from test_dir
# need to make it multiple test cases
path_to_input_file = path_to_test_dir + "test_cases.pkl"

if not os.path.exists(path_to_input_file):
	cut_off_distance = {(1,1):3.7,(1,2):3.7,(2,2):3.7}
	size = 32.130125 - 0.299875
	box_dim = [size, size, size]
	num_of_tests = 2000
	input_param = {'cut_off':cut_off_distance,'box_dim':box_dim,'num_of_tests':num_of_tests}
	print "input_param:", input_param
	with open(path_to_input_file,'w') as f:
		pickle.dump(input_param,f)
else:
	input_param = pickle.load(open(path_to_input_file,'r'))

strain_calculator_run_all_tests(path_to_data_dir, input_param)





	






