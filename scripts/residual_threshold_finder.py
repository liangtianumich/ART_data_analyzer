#!/usr/bin/env python
import os
import multiprocessing as mp
import numpy as np
from correlation_model.correlation_model import residual_threshold_finder

#input
path_to_data_dir = os.environ['DATA_DIR']
# make it a user input by either terminal arguments or from input file
#num_of_tests = 2000
#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [1,2]
num_of_proc = mp.cpu_count()
model = "LinearSVR"
#model = "linear_model"
feature = "displacement"
target =  "shear_strain"
residual_threshold = np.arange(0.01, 1.0, 0.01)
input_param = {"list_of_test_id":list_of_test_id, "num_of_proc": num_of_proc,"model":model,"feature":feature,"target": target,"residual_threshold":residual_threshold}
residual_threshold_finder(path_to_data_dir, input_param)

# results using LinearSVR: at relative threshold:0.54, average number of locally involved atoms: 4.8
