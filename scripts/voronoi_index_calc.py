#!/usr/bin/env python
import os
import pickle
import json
import time
import multiprocessing as mp
from calculator.voronoi_structural_analysis import run_all_tests_voronoi_calculator

path_to_data_dir = os.environ['DATA_DIR']

cut_off = 3.7
l_range = [0.299875, 32.130125]
box_range = [l_range, l_range, l_range]
# num_of_tests = 2000
# list_of_test_id = [x+1 for x in xrange(2000)]
list_of_test_id = [1,2]
num_of_proc = mp.cpu_count()
re_calc = True
atom_list = None
periodic = [True, True, True]
#atom_list = {"local":6}
input_param = {'cut_off':cut_off,'box_range':box_range,'list_of_test_id':list_of_test_id,'num_of_proc':num_of_proc,"re_calc": re_calc, "atom_list":atom_list, "periodic": periodic}

	
if 'num_of_proc' not in input_param:
	input_param["num_of_proc"] = mp.cpu_count()

start_time = time.time()
run_all_tests_voronoi_calculator(path_to_data_dir,input_param)
print "total run time:", time.time() - start_time, "seconds"




	






