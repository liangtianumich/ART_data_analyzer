"""
this event selector module filter each event based on various criterias
each criteria is an independent function
"""
import pandas as pd
import json
import os
import numpy as np
import multiprocessing as mp
from functools import partial
from util import event_energy, Configuration, event_distance, fn_timer



@fn_timer
def event_redudancy_check(path_to_data_dir, input_param, save_results=True, re_calc = False):
	"""
	this function implement stage 2 criteria 3 to remove the redundancy of event pairs
	it will sort all event pair and finally save the unique events into
	a pkl file as a list of strings, such as test1/event_init_sad_fin
	
	criteria 3:
	for the remaining refined searches, any pair is redundant if 
	abs(D(fin - init)_1-D(fin-init)_2) < 0.1 (A)
	AND abs(E(fin-init)_1-E(fin-init))_2 < 0.005(eV)
	AND abs(E(sad-init)_1-E(sad-init))_2 < 0.01(eV)
	"""
	list_of_test_id = input_param["list_of_test_id"]
	
	existed_tests = []
	for i in list_of_test_id:
		path_to_curr_test = path_to_data_dir + "test%s"%i
		if os.path.exists(path_to_curr_test):
			existed_tests.append(i)
	list_of_test_id = existed_tests
	
	box_dim = input_param["box_dim"]
	identical_event_criteria = input_param["identical_event_criteria"]
	if "num_of_proc" in input_param and "re_calc" in input_param:
		num_of_proc = input_param["num_of_proc"]
		re_calc = input_param["re_calc"]
	
	if re_calc is False:
		path_to_final_selected_events = path_to_data_dir + "final_selected_events.json"
		if os.path.exists(path_to_final_selected_events):
			print "reduadancy has already been checked, total number of final selected events:"
			final_events= json.load(open(path_to_final_selected_events,'r'))
			print len(final_events)
			return
			
		
	all_selected_events = get_list_of_selected_events_str(path_to_data_dir, list_of_test_id)
	
	# remove the 2nd identical pair events found
	num_of_selected_events = len(all_selected_events)
	print "total number of selected events:", num_of_selected_events
	removed_index = []
	pool = mp.Pool(processes = num_of_proc)
	for i in xrange(num_of_selected_events):
		if i in removed_index:
			continue
		tests_list = [all_selected_events[j] for j in xrange(i+1,num_of_selected_events)]
		result = pool.map(partial(identical_events,path_to_data_dir=path_to_data_dir,event_1 = all_selected_events[i], box_dim=box_dim,identical_event_criteria=identical_event_criteria), tests_list)
		#for j in xrange(i+1,num_of_selected_events):
		#	is_same = identical_events(path_to_data_dir, all_selected_events[i], all_selected_events[j], box_dim, identical_event_criteria)
		#	if is_same:
		#		removed_index.append(j)
		k=0
		for x in result:
			if x == True:
				removed_index.append(k+i+1)
			k = k+1
	removed_index = np.unique(removed_index)
	# final_selected_events = np.delete(all_selected_events, removed_index).tolist()
	final_selected_events = [i for j, i in enumerate(all_selected_events) if j not in removed_index]
	print "total number of final selected events after removing redundacy:", len(final_selected_events)
	
	#save it into a file called final_selected_events.json
	if save_results is True:
		path_to_final_selected_events = path_to_data_dir + "final_selected_events.json"
		json.dump(final_selected_events, open(path_to_final_selected_events,'w'))
	
	print "done redudancy check for all interested tests!"
	
def get_list_of_selected_events_str(path_to_data_dir, list_of_test_id):
	"""
	this function returns a list containing the strings of the selected events
	after stage I two criteria for all tests in list_of_test_id
	"""
	all_selected_events = []
	for i in list_of_test_id:
		path_to_curr_test = path_to_data_dir + "test%s"%i
		path_to_selected_events = path_to_curr_test + "/results/selected_events.json"
		# skip the test who do not have selected_events.json in all tests specified in
		# list_of_test_id
		if os.path.exists(path_to_selected_events):
			selected_events = json.load(open(path_to_selected_events,'r'))
			# value = [init_state, sad_state, fin_state]
			selected_list = selected_events.values()
			for event in selected_list:
				event_str = ("test%s"%i, [event[0],event[1],event[2]])
				all_selected_events.append(event_str)
	return all_selected_events
@fn_timer
def identical_events(event_2, path_to_data_dir,event_1, box_dim, identical_event_criteria={"D_init_fin": 0.1, "E_init_fin":0.005, "E_init_sad":0.01}):
	"""
	this function return True if two events are identical
	
	criteria 3:
	for the remaining refined searches, any pair is redundant if
	abs(D(fin - init)_1-D(fin-init)_2) < D_init_fin e.g. 0.1 (A) for Cu-Zr
	AND abs(E(fin-init)_1-E(fin-init)_2) < E_init_fin e.g. 0.005(eV) for Cu-Zr
	AND abs(E(sad-init)_1-E(sad-init)_2) < E_init_sad e.g. 0.01(eV) for Cu-Zr
	"""
	D_init_fin = identical_event_criteria["D_init_fin"]
	E_init_fin = identical_event_criteria["E_init_fin"]
	E_init_sad = identical_event_criteria["E_init_sad"]
	
	path_to_event_1_test = path_to_data_dir + event_1[0]
	path_to_event_2_test = path_to_data_dir + event_2[0]
	event_1_init,event_1_sad, event_1_fin = event_1[1][0],event_1[1][1],event_1[1][2]
	event_2_init,event_2_sad, event_2_fin = event_2[1][0],event_2[1][1],event_2[1][2]
	
	event_1_energy = event_energy(path_to_event_1_test)
	
	event_1_init_eng,event_1_sad_eng,event_1_fin_eng = event_1_energy[event_1_init],event_1_energy[event_1_sad],event_1_energy[event_1_fin]
	
	event_2_energy = event_energy(path_to_event_2_test)
	
	event_2_init_eng,event_2_sad_eng,event_2_fin_eng= event_2_energy[event_2_init],event_2_energy[event_2_sad],event_2_energy[event_2_fin]
	
	cond_1 = abs(event_1_fin_eng - event_1_init_eng - (event_2_fin_eng - event_2_init_eng)) < E_init_fin 
	
	cond_2 = abs(event_1_sad_eng - event_1_init_eng - (event_2_sad_eng - event_2_init_eng)) < E_init_sad
	
	if not (cond_1 and cond_2):
		return False
	
	distance_1 = event_distance(path_to_event_1_test, [event_1_init,event_1_sad, event_1_fin],box_dim)
	
	distance_2 = event_distance(path_to_event_2_test, [event_2_init,event_2_sad, event_2_fin],box_dim)
	
	cond_3 = abs(distance_1 -distance_2) < D_init_fin
	
	if cond_3:
		return True
	else:
		return False
	
	#if cond_1 and cond_2 and cond_3:
		#return True
	#else:
		#return False

def df_to_dict(df):
	"""
	this function converts the pandas.Dataframe into a dictionary
	with the key value being the index of pd.Dataframe and the values being
	the values of various columns
	"""
	_final_dict = dict()
	for index,row in df.iterrows():
		_final_dict[index] = row.tolist()
	
	return _final_dict
	
def event_selection(path_to_test_dir = None, box_dim = None, save_results=True, re_calc = False):
	"""
	this function filter all events in a single test, each event has an initial configuration, saddle
	configuration and final configuration. First the event has to be accepted,
	Next this event is selected or not based on whether satisfying the 
	stage 1 two criterias:
	stage 1:
	criteria 1:
	the saddle state energy should be larger than the energy of both initial and final state
	exclude the searches E(sad-init) <0 or E(sad - fin) <0
	
	criteria 2:
	If both the energy difference AND distance between the final state and initial state 
	are small, then the final state is identical to the initial state, should be eliminated
	in details, this happens if dE < 0.02 (eV) AND distance < 1
	
	stage 2 criteria 3 involves each event pair comparison that require much more calculation
	it is implemented as a separate function
	
	criteria 3:
	for the remaining refined searches, any pair is redundant if 
	abs(D(fin - init)_1-D(fin-init)_2) < 0.1 (A)
	AND abs(E(fin-init)_1-E(fin-init))_2 < 0.005(eV)
	AND abs(E(sad-init)_1-E(sad-init))_2 < 0.01(eV)
	Energy stored in log.file.1, with the pattern matching,
	except min1000
	Configuration stored in file min1001
	Total energy Minimum (eV) -8.8897753852E+03
	
	for min1000,
	- Configuration stored in file :           min1000
	__________________________________________________
	- Simulation                   :                 1
	- Attempt                      :                 1
	- Starting from minconf        :              1000
	- Reference Energy (eV)        : -8.8912484758E+03
	
	return:
		selected_events: dict()
			this dict with key being the original event id in event.list
			value being a list containing the string of initial state, saddle state
			final state
	"""
	path_to_selected_events = path_to_test_dir + "/results/selected_events.json"
	if re_calc is False:
		if os.path.exists(path_to_selected_events):
			print "events already selected, load selected_events.json"
			# if selected_event.json is empty, no events selected
			try:
				return json.load(open(path_to_selected_events,'r'))
			except ValueError:
				return None
	print "starting selecting events based on two criteria"
	accepted_events = event_select_accept(path_to_test_dir)
	if accepted_events is None:
		return None
	# event_energy is a function in util.py that extract energy from log.file.1 using regex module
	events_energy = event_energy(path_to_test_dir)
	
	selected_events = dict()
	for x in accepted_events:
		init_state = accepted_events[x][0]
		sad_state = accepted_events[x][1]
		fin_state = accepted_events[x][2]
		list_of_state = [init_state, sad_state, fin_state]
		
		init_energy = events_energy[init_state]
		sad_energy = events_energy[sad_state]
		fin_energy = events_energy[fin_state]
		list_of_energy = [init_energy, sad_energy, fin_energy]
		
		path_to_init_file = path_to_test_dir +'/' + init_state + ".dump"
		path_to_sad_file = path_to_test_dir +'/' + sad_state + ".dump"
		path_to_fin_file = path_to_test_dir +'/' + fin_state + ".dump"
		
		list_of_path = [path_to_init_file, path_to_sad_file, path_to_fin_file]
		
		if single_event_2_criteria(list_of_path, list_of_energy, box_dim):
			selected_events[x] = list_of_state
	
	if selected_events == dict():
		return None
	if save_results is True:
		json.dump(selected_events, open(path_to_selected_events,"w"))
	return selected_events

def event_select_accept(path_to_test_dir = None, save_results=True):
	"""
	this function filter a single event, i.e. initial configuration, saddle
	configuration and final configuration based on whether this event has been accepted
	or not based on the event_list file inside the ART results
	
	this function should only be used in the exe script file, not to be intertwined
	with the calculator module
	
	Input arguments:
		path_to_file: dir
			this path contains the events.list file in current test
			each test has a different path
	return:
		accepted_events: 
			a dictionary with the key being the id of the attempt
			the values being a list containing the initial, saddle, final
			file string
	for example, it will return 
	{0: ['min1000', 'sad1001', 'min1001', 'accepted'], 
	1: ['min1001', 'sad1002', 'min1002', 'accepted'], 
	...
	
	Note: maybe use an extra input argument event_list_df=None as the 2nd 
	option since event_list_df can be passed as both input and output
	to maintain the consistenency of the code for various functions
	"""
	if path_to_test_dir == None:
		raise Exception("no directory path to events.list file has been specified, please specifiy \
		the correct path to the events.list file")
	
	path_to_events_list = path_to_test_dir + "/events.list"
	
	path_to_accepted_events = path_to_test_dir + "/results/accepted_events.json"
	if os.path.exists(path_to_accepted_events):
		return json.load(open(path_to_accepted_events,'r'))
		
	
	events = pd.read_csv(path_to_events_list,sep='\s+', header=None)
	events.columns = ["ini","sad","fin","status"]
	accepted_events = events.loc[events["status"] == "accepted"]
	
	if accepted_events.empty:
		return None
	
	accepted_events = df_to_dict(accepted_events)
	if save_results is True:
		json.dump(accepted_events, open(path_to_accepted_events,"w"))
	# drop the column called "status" since now all events are accepted
	# accepted_events = accepted_events.drop('status', axis=1)
	return accepted_events

def event_3_criteria():
	df = df.loc[df["status"] == "accepted"]
	df["satisfy_3_criteria"] = False
	
	for index,row in df.iterrows():
		if single_event_3_criteria(row) is True:
			df["satisfy_3_criteria"] = True
	return df		

def single_event_2_criteria(list_of_path, list_of_energy, box_dim):
	"""
	this function takes the data of single event, i.e. the configurational data
	and the energy data to check if this event satisfy the stage 1 two criteria
	return True if satisfied, False otherwise
	
	criteria 1:
	the saddle state energy should be larger than the energy of both initial and final state
	exclude the searches E(sad-init) <0 or E(sad - fin) <0
	
	criteria 2:
	If both the energy difference AND distance between the final state and initial state 
	are small, then the final state is identical to the initial state, should be eliminated
	in details, this happens if dE < 0.02 (eV) AND distance < 1
	"""
	
	init_config = Configuration(list_of_path[0],box_dim)
	fin_config = Configuration(list_of_path[2],box_dim)
	
	# criteria 1
	if list_of_energy[1] - list_of_energy[0] < 0 or list_of_energy[1] - list_of_energy[2] < 0:
		return False
	
	#criteria 2
	distance = Configuration.distance_pbc(init_config,fin_config)
	if distance < 1 and abs(list_of_energy[2] - list_of_energy[0]) < 0.02:
		return False
	else:
		return True	
		
	
	
	
	
	
	
