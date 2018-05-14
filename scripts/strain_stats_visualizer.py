#!/usr/bin/env python
import os
from visualizer.strain_visualizer import strain_events_stats_visualization

#input
path_to_data_dir = os.environ['DATA_DIR']
# make it a user input by either terminal arguments or from input file
#num_of_tests = 2000
#list_of_test_id = xrange(num_of_tests+1)
list_of_test_id = [1,2]

# this visualization can also be done by running strain_calc.py with re_calc=False,list_of_test_id = [1,2]
# it takes the events_stats.pkl in each test of list_of_test_id and make the plot
strain_events_stats_visualization(path_to_data_dir, list_of_test_id)
