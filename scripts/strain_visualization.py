#!/usr/bin/env python
import os
import multiprocessing as mp
from visualizer.strain_visualizer import events_strain_visualization

#input
path_to_data_dir = os.environ['DATA_DIR']
# make it a user input by either terminal arguments or from input file
#num_of_tests = 2000
#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [1,2]
num_of_proc = mp.cpu_count()
input_param = {"list_of_test_id":list_of_test_id, "num_of_proc": num_of_proc}
events_strain_visualization(path_to_data_dir, input_param)
