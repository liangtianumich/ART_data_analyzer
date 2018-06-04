#!/usr/bin/env python
import os
import multiprocessing as mp
import numpy as np
from util import run_tests_triggered_atom_is_max_disp

path_to_data_dir = os.environ['DATA_DIR']


#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [1,2]
num_of_proc = mp.cpu_count()

input_param = {"list_of_test_id":list_of_test_id, "num_of_proc": num_of_proc}

run_tests_triggered_atom_is_max_disp(path_to_data_dir, input_param)
