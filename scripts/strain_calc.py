#!/usr/bin/env python
import os
import pickle
import json
import time
import multiprocessing as mp
from calculator.strain_calculator import strain_calculator_run_all_tests_mp

path_to_data_dir = os.environ['DATA_DIR']

# save some backup files that might be useful after running the tests
path_to_test_dir = os.environ['TEST_DIR']

if not os.path.exists(path_to_test_dir):
	os.makedirs(path_to_test_dir)

# make this an input file to load from test_dir
# need to make it multiple test cases
# in the future let strain_calc.py load input file called input.json
path_to_input_file = path_to_test_dir + "/input.json"

cut_off_distance = {(1,1):3.7,(1,2):3.7,(2,2):3.7}
size = 32.130125 - 0.299875
box_dim = [size, size, size]
# num_of_tests = 2000
list_of_test_id = [x+1 for x in xrange(2000)]
num_of_proc = mp.cpu_count()
re_calc = False
atom_list = None
# this will be part of the input file, but not used by strain_calc.py, in 
# event_filter.py
identical_event_criteria = {"D_init_fin": 0.1, "E_init_fin":0.005, "E_init_sad":0.01}
input_param = {'cut_off':cut_off_distance,'box_dim':box_dim,'list_of_test_id':list_of_test_id,'num_of_proc':num_of_proc,"re_calc": re_calc, "atom_list":atom_list,"identical_event_criteria":identical_event_criteria}

	
if 'num_of_proc' not in input_param:
	input_param["num_of_proc"] = mp.cpu_count()

start_time = time.time()
strain_calculator_run_all_tests_mp(path_to_data_dir, input_param)
print "total run time:", time.time() - start_time, "seconds"




	






