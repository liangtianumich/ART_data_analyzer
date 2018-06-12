#!/usr/bin/env python
import os
import multiprocessing as mp
import numpy as np
from correlation_model.correlation_model import all_events_local_atoms_finder

#input
path_to_data_dir = os.environ['DATA_DIR']
# make it a user input by either terminal arguments or from input file
#num_of_tests = 2000
#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [1,2]
num_of_proc = mp.cpu_count()
model = "LinearSVR"
feature = "displacement"
target =  "shear_strain"
re_calc = True
# results using LinearSVR: at relative threshold:0.54, average number of locally involved atoms: 4.8
residual_threshold = 0.54
input_param = {"list_of_test_id":list_of_test_id, "num_of_proc": num_of_proc,"model":model,"feature":feature,"target": target,"re_calc": re_calc}
all_events_local_atoms_finder(path_to_data_dir, input_param, residual_threshold)


