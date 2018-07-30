import numpy as np
import pandas as pd
import os
import multiprocessing as mp
import pickle
import json
from functools import partial
from util import Atom, NN_finder_all, event_strain_disp,event_local_atom_index,read_from_art_input_file, data_dir_to_test_dir
from event_selector import event_selection
from data_reader import *
from visualizer.strain_visualizer import *
from visualizer.general_visualizer import *

def strain_calculator_run_all_tests_mp(path_to_data_dir, input_param):
	"""
	this function run strain calculations on all tests under multiple processors
	"""
	cut_off_distance = input_param["cut_off"]
	box_dim = input_param['box_dim']
	#num_of_tests = input_param['num_of_tests']
	list_of_test_id = input_param['list_of_test_id']
	num_of_proc = input_param['num_of_proc']
	re_calc = input_param["re_calc"]
	# default None calculating all atoms, 
	# if local, make it a dict {"local":8}, only used after calculating all atom strains
	# and perform correlation modeling to find the local number of atoms, then rerun with local mode
	atom_list = input_param["atom_list"]
	
	tests_list = []
	#for i in xrange(num_of_tests+1):
	for i in list_of_test_id:
		try:
			path_to_curr_test = data_dir_to_test_dir(path_to_data_dir, i)
			tests_list.append(path_to_curr_test)
		except Exception:
			pass
	disp_ave, disp_std, disp_max , disp_ave_2, disp_std_2, disp_max_2,disp_ave_3, disp_std_3, disp_max_3 = [], [], [], [], [], [], [], [], []
	
	shear_ave, shear_std, shear_max, shear_ave_2, shear_std_2, shear_max_2,shear_ave_3, shear_std_3, shear_max_3 = [], [], [], [], [], [], [], [], []
	
	vol_ave, vol_std, vol_max, vol_ave_2, vol_std_2, vol_max_2,vol_ave_3, vol_std_3, vol_max_3 = [], [], [], [], [], [], [], [], []
	
	pool = mp.Pool(processes = num_of_proc)
	#partial(my_fun2, general_const=my_const), input_list)
	# search python pool map partial
	# pool.map result preserve the order of input
	result_list = pool.map(partial(strain_calculator_run_single_test,cut_off_distance=cut_off_distance, atom_list=atom_list, box_dim=box_dim, re_calc=re_calc), tests_list)
    
	for curr_test in result_list:
		# if this test do not have data and strain calculation gives None
		# skip this test
		if curr_test is None:
			continue
		for event in curr_test:
			init_sad = event[0]
			sad_fin = event[1]
			init_fin = event[2]
			# calculate the statistics of init_sad and sad_fin		
			disp_ave.append(init_sad["ave"][2])
			disp_std.append(init_sad["std"][2])
			disp_max.append(init_sad["max"][2])
			
			shear_ave.append(init_sad["ave"][1])
			shear_std.append(init_sad["std"][1])
			shear_max.append(init_sad["max"][1])
			
			vol_ave.append(init_sad["ave"][0])
			vol_std.append(init_sad["std"][0])
			vol_max.append(init_sad["max"][0])
			
			disp_ave_2.append(sad_fin["ave"][2])
			disp_std_2.append(sad_fin["std"][2])
			disp_max_2.append(sad_fin["max"][2])
			
			shear_ave_2.append(sad_fin["ave"][1])
			shear_std_2.append(sad_fin["std"][1])
			shear_max_2.append(sad_fin["max"][1])
			
			vol_ave_2.append(sad_fin["ave"][0])
			vol_std_2.append(sad_fin["std"][0])
			vol_max_2.append(sad_fin["max"][0])
			
			disp_ave_3.append(init_fin["ave"][2])
			disp_std_3.append(init_fin["std"][2])
			disp_max_3.append(init_fin["max"][2])
			
			shear_ave_3.append(init_fin["ave"][1])
			shear_std_3.append(init_fin["std"][1])
			shear_max_3.append(init_fin["max"][1])
			
			vol_ave_3.append(init_fin["ave"][0])
			vol_std_3.append(init_fin["std"][0])
			vol_max_3.append(init_fin["max"][0])
			
	pickle.dump({"ave":disp_ave,"std":disp_std,"max":disp_max}, open(path_to_data_dir+"/init_sad_disp_stats.pkl",'w'))
	pickle.dump({"ave":shear_ave,"std":shear_std,"max":shear_max}, open(path_to_data_dir+"/init_sad_shear_stats.pkl",'w'))
	pickle.dump({"ave":vol_ave,"std":vol_std,"max":vol_max}, open(path_to_data_dir+"/init_sad_vol_stats.pkl",'w'))
	
	pickle.dump({"ave":disp_ave_2,"std":disp_std_2,"max":disp_max_2}, open(path_to_data_dir+"/sad_fin_disp_stats.pkl",'w'))
	pickle.dump({"ave":shear_ave_2,"std":shear_std_2,"max":shear_max_2}, open(path_to_data_dir+"/sad_fin_shear_stats.pkl",'w'))
	pickle.dump({"ave":vol_ave_2,"std":vol_std_2,"max":vol_max_2}, open(path_to_data_dir+"/sad_fin_vol_stats.pkl",'w'))
	
	pickle.dump({"ave":disp_ave_3,"std":disp_std_3,"max":disp_max_3}, open(path_to_data_dir+"/init_fin_disp_stats.pkl",'w'))
	pickle.dump({"ave":shear_ave_3,"std":shear_std_3,"max":shear_max_3}, open(path_to_data_dir+"/init_fin_shear_stats.pkl",'w'))
	pickle.dump({"ave":vol_ave_3,"std":vol_std_3,"max":vol_max_3}, open(path_to_data_dir+"/init_fin_vol_stats.pkl",'w'))
	
	
	plot_histogram_3(path_to_data_dir+"/disp_ave.png", [disp_ave,disp_ave_2,disp_ave_3])
	plot_histogram_3(path_to_data_dir+"/disp_std.png", [disp_std,disp_std_2,disp_std_3])
	plot_histogram_3(path_to_data_dir+"/disp_max.png", [disp_max,disp_max_2,disp_max_3])
	
	plot_histogram_3(path_to_data_dir+"/shear_ave.png", [shear_ave,shear_ave_2,shear_ave_3])
	plot_histogram_3(path_to_data_dir+"/shear_std.png", [shear_std,shear_std_2,shear_std_3])
	plot_histogram_3(path_to_data_dir+"/shear_max.png", [shear_max,shear_max_2,shear_max_3])
	
	plot_histogram_3(path_to_data_dir+"/vol_ave.png", [vol_ave,vol_ave_2,vol_ave_3])
	plot_histogram_3(path_to_data_dir+"/vol_std.png", [vol_std,vol_std_2,vol_std_3])
	plot_histogram_3(path_to_data_dir+"/vol_max.png", [vol_max,vol_max_2,vol_max_3])	
	print "done!"

def strain_calculator_run_single_test(test, cut_off_distance, box_dim, atom_list = None, save_results = True, re_calc=False):
	"""
	this function perform strain and displacement calculation for a single test
	whose dir is test
	
	this atom_list would be the item_id of local atoms around the triggered atom,
	which should be event dependent. This means that after we get the number of involved
	local atoms, we need to write a function to extract the local atoms index around
	the triggered atom, which will be our atom_list for this event for local event,
	however, this will not be the issue for the test containing single event
	"""
	path_to_curr_result = test + "/results"
	if not os.path.exists(path_to_curr_result):
		os.makedirs(path_to_curr_result)

	# get each of the selected events for current test
	path_to_event_list = path_to_curr_result + "/selected_events.json"
	#escape when this test do not have log.file.1 file or configuration file
	try:
		event_list = event_selection(test,box_dim,re_calc = re_calc)
	except IOError:
		return None
	# if no events has been selected in this test
	if event_list is None:
		return None
	
	# if final_selected_events.json exists, this means that event redudancy check
	# happens before calculating all strains
	# check if each event is in the final_selected_events.json file
	path_to_data_dir,test_id = os.path.split(test)
	if 'test' not in test_id:
		test_id = 'test%s'%test_id
	
	path_to_final_events = path_to_data_dir + "/final_selected_events.json"
	event_list = event_list.values()
	
	
	if os.path.exists(path_to_final_events):
		target_events = []
		final_events = json.load(open(path_to_final_events,'r'))
		for event in final_events:
			if event[0] == test_id:
				if event[1] in event_list:
					target_events.append(event[1])
		event_list = target_events
					
	
	# list in the order of vol_strain, shear_strain, disp
	test_results = []
	
	test_atom_list = atom_list
	# for each event, init to sad and sad to fin
	for event in event_list:
		init_sad_event_result = dict()
		sad_fin_event_result = dict()
		init_fin_event_result = dict()
		
		init, sad, fin = event[0], event[1], event[2]
		path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
		print ('\n')
		print "path_to_curr_event", path_to_curr_event
		print ('\n')
		if not os.path.exists(path_to_curr_event):
			os.makedirs(path_to_curr_event)
		path_to_file_ini = test + '/' + init + ".dump"
		path_to_file_sad = test + '/' + sad + ".dump"
		path_to_file_fin = test + '/' + fin + ".dump"
		
		initial_config_data = read_data_from_file(path_to_file_ini)
		saddle_config_data = read_data_from_file(path_to_file_sad)
		final_config_data = read_data_from_file(path_to_file_fin)
		
		path_to_init_sad = path_to_curr_event + "/init_sad"
		path_to_sad_fin = path_to_curr_event + "/sad_fin"
		path_to_init_fin = path_to_curr_event + "/init_fin"
		
		if not os.path.exists(path_to_init_sad):
			os.makedirs(path_to_init_sad)
		if not os.path.exists(path_to_sad_fin):
			os.makedirs(path_to_sad_fin)
		if not os.path.exists(path_to_init_fin):
			os.makedirs(path_to_init_fin)		
		
		if test_atom_list is None:
			test_atom_list = (initial_config_data["item"]).tolist()
		
		# if atom_list is a dict,e.g. {"local":4}
		# find the atom list of item_id of 4 local atoms around triggered atom
		if type(test_atom_list) == dict or test_atom_list == "local":
			print "\n starting local mode strain calculations"
			if type(test_atom_list) == dict:
				triggered_atom_index = read_from_art_input_file(test)
				num_of_involved_atom = test_atom_list["local"]
				# triggered atom is only related to initial config, for sad_fin, 
				# also use the initial config to find the local NN index
				atom_list = event_local_atom_index(initial_config_data, triggered_atom_index, num_of_involved_atom, path_to_curr_event, box_dim, re_calc=re_calc)
				print "NN item_id list (distance in increasing order including triggering atoms) of the triggering atoms"
			elif test_atom_list == "local":
				path_to_local_atom_index = path_to_curr_event + "/local_atoms_index.json"
				if os.path.exists(path_to_local_atom_index):
					local_atom_list = json.load(open(path_to_local_atom_index,"r"))
					atom_list = [atom + 1 for atom in local_atom_list]
				else:
					raise Exception("local_atoms_index.json file not exists, run local_atoms_index_finder.py first")
			
			print "in initial configuration that will be used for all local strain calculations:", atom_list
			if atom_list == []:
				continue
			print "\n initial to saddle: \n"
			init_sad_strain,init_sad_disp = local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_init_sad, atom_list = atom_list,local=True,re_calc = re_calc)
			print "\n saddle to final: \n"
			sad_fin_strain,sad_fin_disp = local_strain_calculator_orth(saddle_config_data, final_config_data, cut_off_distance, box_dim, path_to_sad_fin, atom_list = atom_list,local=True,re_calc = re_calc)
			print "\n initial to final: \n"
			init_fin_strain,init_fin_disp = local_strain_calculator_orth(initial_config_data, final_config_data, cut_off_distance, box_dim, path_to_init_fin, atom_list = atom_list,local=True,re_calc = re_calc)
		else:
			print "\n initial to saddle: \n"
			init_sad_strain,init_sad_disp = local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_init_sad, atom_list = atom_list, re_calc = re_calc)
			print "\n saddle to final: \n"
			sad_fin_strain,sad_fin_disp = local_strain_calculator_orth(saddle_config_data, final_config_data, cut_off_distance, box_dim, path_to_sad_fin, atom_list = atom_list, re_calc = re_calc)
			print "\n initial to final: \n"
			init_fin_strain,init_fin_disp = local_strain_calculator_orth(initial_config_data, final_config_data, cut_off_distance, box_dim, path_to_init_fin, atom_list = atom_list, re_calc = re_calc)
			
		#init_sad_vol_strain, init_sad_shear_strain, init_sad_displacement = event_strain_disp(init_sad_strain,init_sad_disp)
		#sad_fin_vol_strain, sad_fin_shear_strain, sad_fin_displacement = event_strain_disp(sad_fin_strain,sad_fin_disp)
		init_sad = event_strain_disp(init_sad_strain,init_sad_disp)
		sad_fin = event_strain_disp(sad_fin_strain,sad_fin_disp)
		init_fin = event_strain_disp(init_fin_strain,init_fin_disp)
		
		
		init_sad_event_result['ave']=[np.mean(init_sad[0]),np.mean(init_sad[1]),np.mean(init_sad[2])]
		init_sad_event_result['std']=[np.std(init_sad[0]),np.std(init_sad[1]),np.std(init_sad[2])]
		init_sad_event_result['max']=[np.max(init_sad[0]),np.max(init_sad[1]),np.max(init_sad[2])]
		
		sad_fin_event_result['ave']=[np.mean(sad_fin[0]),np.mean(sad_fin[1]),np.mean(sad_fin[2])]
		sad_fin_event_result['std']=[np.std(sad_fin[0]),np.std(sad_fin[1]),np.std(sad_fin[2])]
		sad_fin_event_result['max']=[np.max(sad_fin[0]),np.max(sad_fin[1]),np.max(sad_fin[2])]
		
		init_fin_event_result['ave']=[np.mean(init_fin[0]),np.mean(init_fin[1]),np.mean(init_fin[2])]
		init_fin_event_result['std']=[np.std(init_fin[0]),np.std(init_fin[1]),np.std(init_fin[2])]
		init_fin_event_result['max']=[np.max(init_fin[0]),np.max(init_fin[1]),np.max(init_fin[2])]
		
		event_state = [test_id,[init,sad,fin]]
		
		test_results.append((init_sad_event_result,sad_fin_event_result,init_fin_event_result,event_state))
	
	if save_results is True:
		print "begin saving current test each event statistics results into pkl file"
		with open(path_to_curr_result +"/events_stats.pkl", 'w') as f:
			pickle.dump(test_results,f)
			f.close()
	return test_results

def local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_test_dir, atom_list = None, local=False, save_results = True, re_calc = False):
	"""
	this function calculate various local atomic strain quantities for atoms whose item id stored inside
	the atom_list under periodic boundary condition for orthogonal simulation box
	input arguments:
	
	initial_config_data: instance of pandas.Dataframe
		intial configuration data
	
	saddle_config_data: instance of pandas.Dataframe
		saddle configuration data
		
	cut_off_distance: dict
		dictionary contains the multiple pair cut-off distance
		currently use tuples as keys for immutability
		for example, {(1,1):3.7,(1,2):2.7,(2,2):3.0} means that atom_type 1 and 1 cut-off
		is 3.7, etc
		It is suggested that cut_off_distance are used as a user input,
		not implemented by calculating from the first min of pair distribution function
		of initial_config_data since it is not a module that reused a lot
		it is good enough to get this single value through Ovide coordination analysis rdf or other software
		instead of a rigourous calculation
	
	box_dim: list,
		the spatial dimension of simulation box in [x,y,z]
	
	path_to_test_dir: str
		str of directory path to the test results dir, under it, it will save the data into strain_results.pkl
	
	atom_list: list,
		a list storing the atom item id of interested atoms
		default calculating all atoms in the initial_config_data
	
	save_results: boolean, default True
		if True, save calculated local strains results into a pkl file
		
	returns:
		strain: dict()
			a dictionary with key being the item id number of the interested atom
			values being the calculated various local atomic strains measures 
			stored inside in a list
		disp_results: dict()
			a dictionary with key being the item id number of the itnerested atoms
			value being the calcualted displacement for this atom under pbc
			
	Note:
		the nearest neighbor is determined through the atomic configuration in
		initial_config_data
	"""
	
	# check if the strain_results_dict.pkl nn_results_dict.pkl file exists or not
	if local is False:
		path_to_strain_results = path_to_test_dir + "/strain_results_dict.pkl"
		path_to_displacement = path_to_test_dir + "/displacement_results_dict.pkl"
	elif local is True:
		path_to_all_strain_results = path_to_test_dir + "/strain_results_dict.pkl"
		path_to_all_displacement = path_to_test_dir + "/displacement_results_dict.pkl"
		
		path_to_strain_results = path_to_test_dir + "/local_strain_results_dict.pkl"
		path_to_displacement = path_to_test_dir + "/local_displacement_results_dict.pkl"
	else:
		raise Exception("local mode is either True or False")
	
	if re_calc is False:
		if os.path.exists(path_to_strain_results) and os.path.exists(path_to_displacement):
			print "atomic strain and displacement has already been calculated and saved in pkl file, skip"
			return (pickle.load(open(path_to_strain_results,'r')), pickle.load(open(path_to_displacement,'r')))
		if local is True:
			if os.path.exists(path_to_all_strain_results) and os.path.exists(path_to_all_displacement):
				all_strains = pickle.load(open(path_to_all_strain_results,'r'))
				all_disp = pickle.load(open(path_to_all_displacement,'r'))
				local_strains = dict((k, all_strains[k]) for k in atom_list if k in all_strains)
				local_disp = dict((k, all_disp[k]) for k in atom_list if k in all_disp)
				return (local_strains, local_disp)

	print "starting calculating atomic strain and displacement magnitude"
	strain = dict()
	disp_results = dict()

	# if not specifying the atom_list, choose all atoms in the initial_config_data
	if atom_list is None:
		atom_list = (initial_config_data["item"]).tolist()
	_data = initial_config_data
	
	nn = NN_finder_all(initial_config_data, cut_off_distance, box_dim, path_to_test_dir, atom_list, re_calc = re_calc)
	
	for item in atom_list:
		#calcualte displacement
		init_atom = from_ds(initial_config_data.loc[initial_config_data["item"]==item])
		#init_atom.box_dim = box_dim
		sad_atom = from_ds(saddle_config_data.loc[saddle_config_data["item"]==item])
		#sad_atom.box_dim = box_dim
		disp_results[item] = distance_pbc(init_atom, sad_atom, box_dim)
		
		# calcualte strain
		item_nn = nn[item]
		NN_list = item_nn["item"].tolist()
		NN_list.append(item)
		# NN_intial is the pandas.Dataframe of intial data with interested_atom and its NN
		NN_initial = initial_config_data.loc[initial_config_data["item"].isin(NN_list)]
		
		# NN_saddle is the pandas.Dataframe of saddle data with interested_atom and its NN
		NN_saddle = saddle_config_data.loc[saddle_config_data["item"].isin(NN_list)]
		#NN_intial = initial_config_data.loc[(initial_config_data['item'].isin(item_nn["item"].tolist()))]
		#NN_intial.append(item)
		#NN_saddle = (saddle_config_data['item'].isin(item_nn["item"])).tolist() 
		#saddle_config_nn = saddle_config_data.loc[(saddle_config_data['item'].isin(item_nn["item"])).tolist()]
		
		# local_strains should be a dict as well since it has multiple output strain
		# or a list
		local_strains = local_strain_calculator_atom_orth(NN_initial, NN_saddle, item, box_dim)
		strain[item] = local_strains
			
	
	if save_results is True:
		print "begin saving atomic strain and displacement dict results into pkl file"
		with open(path_to_strain_results, 'w') as f:
			pickle.dump(strain,f)
			f.close()
		with open(path_to_displacement, 'w') as f:
			pickle.dump(disp_results,f)
			f.close()
		print "atomic strain and displacement results saved into pkl file as dictionary"
	print "done strain and displacement calculations!"
	return (strain, disp_results)


def local_strain_calculator_atom_orth(initial_config_atom, saddle_config_atom, atom_item, box_dim):
	"""
	this function calculate the local atomic strain according to the algorithm
	proposed in http://li.mit.edu/Archive//Graphics/A/annotate_atomic_strain/Doc/main.pdf
	
	Input arguments:
		
		initial_config_atom: an instance of pandas.Dataframe
			interested atom and its NN atoms in intial configuration data stores atom id, atom type and x,y,z coordinates of the dump file
		
		saddle_config_atom: an instance of pandas.Dataframe
			interested atom and its NN atoms in saddle configuration data stores atom id, atom type and x,y,z coordinates of the dump file
		
		atom_item: integer
			the item id of interested atom
			
		box_dim: list
			the simulation box spatial dimension to implement periodic boundary condition
	return:
		atom_strain: list
			storing various local atomic strains in list
			such as the von mises shear strain invariants, hydrostatic invariant 
	"""
	# calculate dij and d0_ij
	NN_initial = initial_config_atom.loc[initial_config_atom['item'] != atom_item]
	Atom_initial = initial_config_atom.loc[initial_config_atom['item'] == atom_item]
	Atom_ini_obj = from_ds(Atom_initial)
	
	
	NN_saddle = saddle_config_atom.loc[saddle_config_atom['item'] != atom_item]
	Atom_saddle = saddle_config_atom.loc[saddle_config_atom['item'] == atom_item]
	Atom_sad_obj = from_ds(Atom_saddle)
	
	Dim = 3
	V = np.zeros(shape=(Dim,Dim))
	W = np.zeros(shape=(Dim,Dim))
	
	for (index,atom_ini_NN) in NN_initial.iterrows():
		# d0_ji in pandas.Series
		d0_ji = from_ds(atom_ini_NN) - Atom_ini_obj
		#d0_ji = to_list(d0_ji)
		atom_sad_NN = from_ds(NN_saddle.loc[NN_saddle["item"] == atom_ini_NN["item"]])
		# d_ji in pandas.Series
		d_ji = atom_sad_NN - Atom_sad_obj
		#d_ji = to_list(d_ji)
		
		# begin implement pbc for d_ji and d0_ji as in 
		# https://en.wikipedia.org/wiki/Periodic_boundary_conditions
		# dx = x[j] - x[i];
		# dx -= x_size * nearbyint(dx * x_rsize)
		
		Near_int_d = [int(round(i)) for i in np.array(d_ji) * 1.0/np.array(box_dim)]
		d_ji = np.array(d_ji) - np.array(box_dim) * np.array(Near_int_d)
		
		Near_int_d0 = [int(round(i)) for i in np.array(d0_ji) * 1.0/np.array(box_dim)]
		d0_ji = np.array(d0_ji) - np.array(box_dim) * Near_int_d0
		
		
		# begin calculate the V and M matrix
		V[0][0] = V[0][0] + d0_ji[0] * d0_ji[0]
		V[0][1] = V[0][1] + d0_ji[0] * d0_ji[1]
		V[0][2] = V[0][2] + d0_ji[0] * d0_ji[2]
		V[1][0] = V[1][0] + d0_ji[1] * d0_ji[0]
		V[1][1] = V[1][1] + d0_ji[1] * d0_ji[1]
		V[1][2] = V[1][2] + d0_ji[1] * d0_ji[2]
		V[2][0] = V[2][0] + d0_ji[2] * d0_ji[0]
		V[2][1] = V[2][1] + d0_ji[2] * d0_ji[1]
		V[2][2] = V[2][2] + d0_ji[2] * d0_ji[2]

		W[0][0] = W[0][0] + d0_ji[0] * d_ji[0]
		W[0][1] = W[0][1] + d0_ji[0] * d_ji[1]
		W[0][2] = W[0][2] + d0_ji[0] * d_ji[2]
		W[1][0] = W[1][0] + d0_ji[1] * d_ji[0]
		W[1][1] = W[1][1] + d0_ji[1] * d_ji[1]
		W[1][2] = W[1][2] + d0_ji[1] * d_ji[2]
		W[2][0] = W[2][0] + d0_ji[2] * d_ji[0]
		W[2][1] = W[2][1] + d0_ji[2] * d_ji[1]
		W[2][2] = W[2][2] + d0_ji[2] * d_ji[2]
	
	J = np.dot(np.linalg.inv(V), W)
	
	mu = (np.dot(J, np.transpose(J)) - np.identity(3)) * 0.5
	
	mu_hydro = np.trace(mu)/3.0
	
	mu_off_diag = (mu - mu_hydro * np.identity(3))
	
	mu_Mises = (0.5 * np.trace(np.dot(mu_off_diag, mu_off_diag))) ** 0.5
	
	return [mu_hydro, mu_Mises]

def from_ds(data):
	"""
	used inside iteration for run time reduction
	create a Atom class object from single column pandas.DataFrame or pandas.Series
	Input argument:
		data: instance of single column pandas.DataFrame or pandas.Series
			data of a single atom
	
	return:
		a numpy array of atom location
	"""
	_data = data.squeeze()
	return np.array((_data["x"], _data["y"], _data["z"]))

def distance_pbc(atom_1, atom_2, box_dim):
	"""
	use inside large iterations to improve code performance
	method to calculate the distance between atom_1 and atom_2 under
	periodic boundary condition, where atom_1 and atom_2 are the 
	instances of class Atoms
	
	Arguments:
		atom_1: instance of class Atom 
			atomic coordinates of atom_1
		atom_2: instance of class Atom
			atomic coordinates of atom_2
	return:
		distance: float
			the minimum pair distance between atom_1 and atom_2 under pbc
	"""
	lx,ly,lz = box_dim[0],box_dim[1],box_dim[2]
	
	_pair_list = np.array(([0,0,0],[lx,0,0],[-lx,0,0],[0,ly,0],[0,-ly,0],[0,0,lz],[0,0,-lz]))
	
	return min([np.linalg.norm(atom_1 - atom_2 + _pair) for _pair in _pair_list])


