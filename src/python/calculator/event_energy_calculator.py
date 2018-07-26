"""
this event energy calculator module extract various energy related quantities from
each event
"""
from util import event_act_relax_energy, get_list_of_selected_events_str
import multiprocessing as mp
import pandas as pd
import numpy as np
import os, json
from visualizer.event_energy_visualizer import plot_act_relax_histogram

def energy_calculator_run_all_tests_mp(path_to_data_dir, input_param, save_results = True, re_calc = False):
	"""
	this function extract the activation energy and relaxation energy
	for all final filtered events in all tests in list_of_test_id, summary of their statistics
	, for the purpose of checking the convergence of their distribution
	"""
	list_of_test_id = input_param['list_of_test_id']
	num_of_proc = input_param['num_of_proc']
	re_calc = input_param["re_calc"]
	path_to_all_act_relax_eng = path_to_data_dir + "act_relax_eng_filtered_events.json"
	
	if re_calc is False:
		if os.path.exists(path_to_all_act_relax_eng):
			saved_results = json.load(open(path_to_all_act_relax_eng, 'r'))
			all_act_energy, all_relax_eng = [], []
			for result in saved_results:
				all_act_energy.append(result[1])
				all_relax_eng.append(result[2])
			return [all_act_eng, all_relax_eng]
			
		
	
	list_event_str = get_list_of_final_filtered_events_str(path_to_data_dir)
	pool = mp.Pool(processes = num_of_proc)
	result = pool.map(partial(event_act_relax_energy,path_to_data_dir = path_to_data_dir), list_event_str)
	
	all_act_eng, all_relax_eng = [], []
	saved_data = []
	for event_result in result:
		(event, act_eng, relax_eng) = event_result
		saved_data.append([event,act_eng,relax_eng])
		all_act_eng.append(act_eng)
		all_relax_eng.append(relax_eng)
	
	# get statistics
	energy_results = pd.DataFrame({"act_eng":np.array(all_act_eng),"relax_eng": np.array(all_relax_eng)})
	energy_results.describe()
	
	# plot the energy distribution
	path_to_eng_plot = path_to_data_dir + "act_relax_eng_filtered_events.png"
	
	plot_act_relax_histogram(path_to_eng_plot, [all_act_eng, all_relax_eng])
	print "done finding all activation energy and relaxation energy for all filtered events in list_of_test_id"
	
	if save_results is True:
		json.dump(saved_data, open(path_to_all_act_relax_eng,'w'))

	return [all_act_eng, all_relax_eng]
	

	
	
		
		
		
