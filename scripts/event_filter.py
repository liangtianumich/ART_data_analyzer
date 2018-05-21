#!/usr/bin/env python
import os
from event_selector import event_redudancy_check
import multiprocessing as mp

#input
path_to_data_dir = os.environ['DATA_DIR']
# make it a user input by either terminal arguments or from input file
#num_of_tests = 2000
#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [x+1 for x in xrange(2000)]

size = 32.130125 - 0.299875
box_dim = [size, size, size]
identical_event_criteria = {"D_init_fin": 0.1, "E_init_fin": 0.005, "E_init_sad": 0.01}

num_of_proc = mp.cpu_count()

re_calc = True

# identical_event_criteria
input_param = {"list_of_test_id": list_of_test_id, "box_dim": box_dim, "identical_event_criteria": identical_event_criteria,"num_of_proc": num_of_proc, "re_calc": re_calc}

event_redudancy_check(path_to_data_dir, input_param)
	
	
