import numpy as np
import pandas as pd
import os
import pickle
import json
from util import Atom, NN_finder_all
from event_selector import event_selection
from data_reader import *
from visualizer.strain_visualizer import *

def event_strain_disp(event_strain_dict,event_disp_dict):
	"""
	this function takes a strain dictionary with the key being atom item_id
	and value being a list of [volume, von_Mises] converts it into a series statistics
	"""
	vol_strain = []
	shear_strain = []
	disp = []
	
	keys = list(event_strain_dict.keys())
	keys.sort()
	for i in keys:
		vol_strain.append(event_strain_dict[i][0])
		shear_strain.append(event_strain_dict[i][1])
		disp.append(event_disp_dict[i])
	return (vol_strain, shear_strain, disp)

def strain_calculator_run_all_tests(path_to_data_dir, input_param, re_calc = False):
	"""
	this function run all tests starting with test* inside a data directory
	path_to_input = path_to_curr_result + "/input.json"
	"""
	cut_off_distance = input_param["cut_off"]
	box_dim = input_param['box_dim']
	num_of_tests = input_param['num_of_tests']
	
	tests_list = []
	for i in xrange(num_of_tests+1):
		path_to_curr_test = path_to_data_dir + "test%s"%i
		if os.path.exists(path_to_curr_test):
			tests_list.append(path_to_curr_test)
	
	disp_ave, disp_std, disp_max , disp_ave_2, disp_std_2, disp_max_2 = [], [], [], [], [], []
	
	shear_ave, shear_std, shear_max, shear_ave_2, shear_std_2, shear_max_2 = [], [], [], [], [], []
	
	vol_ave, vol_std, vol_max, vol_ave_2, vol_std_2, vol_max_2 = [], [], [], [], [], []
	
	
	for test in tests_list:
		# path to current test results dir
		path_to_curr_result = test + "/results"
		if not os.path.exists(path_to_curr_result):
			os.makedirs(path_to_curr_result)

		# get each of the selected events for current test
		path_to_event_list = path_to_curr_result + "/selected_events.json"
		
		event_list = event_selection(test,box_dim,re_calc = re_calc)
		
		# for each event, init to sad and sad to fin
		for (index,event) in event_list.items():		
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
			
			if not os.path.exists(path_to_init_sad):
				os.makedirs(path_to_init_sad)
			if not os.path.exists(path_to_sad_fin):
				os.makedirs(path_to_sad_fin)
			
			print "\n initial to saddle: \n"
			init_sad_strain,init_sad_disp = local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_init_sad, re_calc = re_calc)
			print "\n saddle to final: \n"
			sad_fin_strain,sad_fin_disp = local_strain_calculator_orth(saddle_config_data, final_config_data, cut_off_distance, box_dim, path_to_sad_fin, re_calc = re_calc)
			
			#init_sad_vol_strain, init_sad_shear_strain, init_sad_displacement = event_strain_disp(init_sad_strain,init_sad_disp)
			#sad_fin_vol_strain, sad_fin_shear_strain, sad_fin_displacement = event_strain_disp(sad_fin_strain,sad_fin_disp)
			init_sad = event_strain_disp(init_sad_strain,init_sad_disp)
			sad_fin = event_strain_disp(sad_fin_strain,sad_fin_disp)
			
			path_to_init_sad_disp_strain = path_to_init_sad + '/disp_shear_strain.png'
			plot_2d_shear(path_to_init_sad_disp_strain,init_sad[2], init_sad[1])
			
			path_to_init_sad_disp_vol_strain = path_to_init_sad + '/disp_vol_strain.png'
			plot_2d_vol(path_to_init_sad_disp_vol_strain,init_sad[2], init_sad[0])
			
			path_to_sad_fin_disp_strain = path_to_sad_fin + '/disp_shear_strain.png'
			plot_2d_shear(path_to_sad_fin_disp_strain,sad_fin[2], sad_fin[1])
			
			path_to_sad_fin_disp_vol_strain = path_to_sad_fin + '/disp_vol_strain.png'
			plot_2d_vol(path_to_sad_fin_disp_vol_strain,init_sad[2], init_sad[0])
			
			pickle.dump(init_sad[2], open(path_to_init_sad+"/disp_results_list.pkl",'w'))
			pickle.dump(init_sad[1], open(path_to_init_sad+"/shear_strain_results_list.pkl",'w'))
			pickle.dump(init_sad[0], open(path_to_init_sad+"/vol_strain_results_list.pkl",'w'))
			
			pickle.dump(sad_fin[2], open(path_to_sad_fin+"/disp_results_list.pkl",'w'))
			pickle.dump(sad_fin[1], open(path_to_sad_fin+"/shear_strain_results_list.pkl",'w'))
			pickle.dump(sad_fin[0], open(path_to_sad_fin+"/vol_strain_results_list.pkl",'w'))
			
			plot_histogram_2(path_to_curr_event + "/disp_histogram.png", [init_sad[2],sad_fin[2]])
			plot_histogram_2(path_to_curr_event + "/shear_strain_histogram.png", [init_sad[1],sad_fin[1]])
			plot_histogram_2(path_to_curr_event + "/vol_strain_histogram.png", [init_sad[0],sad_fin[0]])
			
			#plot_histogram(path_to_init_sad + "/disp_histogram.tif", init_sad[2])
			#plot_histogram(path_to_init_sad + "/shear_strain_histogram.tif", init_sad[1])
			#plot_histogram(path_to_init_sad + "/vol_strain_histogram.tif", init_sad[0])
			
			#plot_histogram(path_to_sad_fin + "/disp_histogram.tif", sad_fin[2])
			#plot_histogram(path_to_sad_fin + "/shear_strain_histogram.tif", sad_fin[1])
			#plot_histogram(path_to_sad_fin + "/vol_strain_histogram.tif", sad_fin[0])
			
			# calculate the statistics of init_sad and sad_fin		
			disp_ave.append(np.mean(init_sad[2]))
			disp_std.append(np.std(init_sad[2]))
			disp_max.append(np.max(init_sad[2]))
			
			shear_ave.append(np.mean(init_sad[1]))
			shear_std.append(np.std(init_sad[1]))
			shear_max.append(np.max(init_sad[1]))
			
			vol_ave.append(np.mean(init_sad[0]))
			vol_std.append(np.std(init_sad[0]))
			vol_max.append(np.max(init_sad[0]))
			
			
			disp_ave_2.append(np.mean(sad_fin[2]))
			disp_std_2.append(np.std(sad_fin[2]))
			disp_max_2.append(np.max(sad_fin[2]))
			
			shear_ave_2.append(np.mean(sad_fin[1]))
			shear_std_2.append(np.std(sad_fin[1]))
			shear_max_2.append(np.max(sad_fin[1]))
			
			vol_ave_2.append(np.mean(sad_fin[0]))
			vol_std_2.append(np.std(sad_fin[0]))
			vol_max_2.append(np.max(sad_fin[0]))
	
	
	pickle.dump({"ave":disp_ave,"std":disp_std,"max":disp_max}, open(path_to_data_dir+"/init_sad_disp_stats.pkl",'w'))
	pickle.dump({"ave":shear_ave,"std":shear_std,"max":shear_max}, open(path_to_data_dir+"/init_sad_shear_stats.pkl",'w'))
	pickle.dump({"ave":vol_ave,"std":vol_std,"max":vol_max}, open(path_to_data_dir+"/init_sad_vol_stats.pkl",'w'))
	
	pickle.dump({"ave":disp_ave_2,"std":disp_std_2,"max":disp_max_2}, open(path_to_data_dir+"/sad_fin_disp_stats.pkl",'w'))
	pickle.dump({"ave":shear_ave_2,"std":shear_std_2,"max":shear_max_2}, open(path_to_data_dir+"/sad_fin_shear_stats.pkl",'w'))
	pickle.dump({"ave":vol_ave_2,"std":vol_std_2,"max":vol_max_2}, open(path_to_data_dir+"/sad_fin_vol_stats.pkl",'w'))
	
	#pickle.dump(disp_ave, open(path_to_data_dir+"/init_sad_disp_ave.pkl",'w'))
	#pickle.dump(disp_std, open(path_to_data_dir+"/init_sad_disp_std.pkl",'w'))
	#pickle.dump(disp_max, open(path_to_data_dir+"/init_sad_disp_max.pkl",'w'))
	
	#pickle.dump(shear_ave, open(path_to_data_dir+"/init_sad_shear_ave.pkl",'w'))
	#pickle.dump(shear_std, open(path_to_data_dir+"/init_sad_shear_std.pkl",'w'))
	#pickle.dump(shear_max, open(path_to_data_dir+"/init_sad_shear_max.pkl",'w'))
	
	#pickle.dump(vol_ave, open(path_to_data_dir+"/init_sad_vol_ave.pkl",'w'))
	#pickle.dump(vol_std, open(path_to_data_dir+"/init_sad_vol_std.pkl",'w'))
	#pickle.dump(vol_max, open(path_to_data_dir+"/init_sad_vol_max.pkl",'w'))
	
	#pickle.dump(disp_ave_2, open(path_to_data_dir+"/sad_init_disp_ave.pkl",'w'))
	#pickle.dump(disp_std_2, open(path_to_data_dir+"/sad_init_disp_std.pkl",'w'))
	#pickle.dump(disp_max_2, open(path_to_data_dir+"/sad_init_disp_max.pkl",'w'))
	
	#pickle.dump(shear_ave_2, open(path_to_data_dir+"/sad_init_shear_ave.pkl",'w'))
	#pickle.dump(shear_std_2, open(path_to_data_dir+"/sad_init_shear_std.pkl",'w'))
	#pickle.dump(shear_max_2, open(path_to_data_dir+"/sad_init_shear_max.pkl",'w'))
	
	#pickle.dump(vol_ave_2, open(path_to_data_dir+"/sad_init_vol_ave.pkl",'w'))
	#pickle.dump(vol_std_2, open(path_to_data_dir+"/sad_init_vol_std.pkl",'w'))
	#pickle.dump(vol_max_2, open(path_to_data_dir+"/sad_init_vol_max.pkl",'w'))
	
	plot_histogram_2(path_to_data_dir+"/disp_ave.png", [disp_ave,disp_ave_2])
	plot_histogram_2(path_to_data_dir+"/disp_std.png", [disp_std,disp_std_2])
	plot_histogram_2(path_to_data_dir+"/disp_max.png", [disp_max,disp_max_2])
	
	plot_histogram_2(path_to_data_dir+"/shear_ave.png", [shear_ave,shear_ave_2])
	plot_histogram_2(path_to_data_dir+"/shear_std.png", [shear_std,shear_std_2])
	plot_histogram_2(path_to_data_dir+"/shear_max.png", [shear_max,shear_max_2])
	
	plot_histogram_2(path_to_data_dir+"/vol_ave.png", [vol_ave,vol_ave_2])
	plot_histogram_2(path_to_data_dir+"/vol_std.png", [vol_std,vol_std_2])
	plot_histogram_2(path_to_data_dir+"/vol_max.png", [vol_max,vol_max_2])
	
			
	#plot_histogram(path_to_data_dir+"/init_sad_disp_ave.tif", disp_ave)
	#plot_histogram(path_to_data_dir+"/init_sad_disp_std.tif", disp_std)
	#plot_histogram(path_to_data_dir+"/init_sad_disp_max.tif", disp_max)
	
	#plot_histogram(path_to_data_dir+"/sad_fin_disp_ave.tif", disp_ave_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_disp_std.tif", disp_std_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_disp_max.tif", disp_max_2)
	
	#plot_histogram(path_to_data_dir+"/init_sad_shear_ave.tif", shear_ave)
	#plot_histogram(path_to_data_dir+"/init_sad_shear_std.tif", shear_std)
	#plot_histogram(path_to_data_dir+"/init_sad_shear_max.tif", shear_max)
	
	#plot_histogram(path_to_data_dir+"/sad_fin_shear_ave.tif", shear_ave_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_shear_std.tif", shear_std_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_shear_max.tif", shear_max_2)
	
	#plot_histogram(path_to_data_dir+"/init_sad_vol_ave.tif", vol_ave)
	#plot_histogram(path_to_data_dir+"/init_sad_vol_std.tif", vol_std)
	#plot_histogram(path_to_data_dir+"/init_sad_vol_max.tif", vol_max)
	
	#plot_histogram(path_to_data_dir+"/sad_fin_vol_ave.tif", vol_ave_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_vol_std.tif", vol_std_2)
	#plot_histogram(path_to_data_dir+"/sad_fin_vol_max.tif", vol_max_2)
	
	print "done!"

def local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_test_dir, atom_list = None, save_results = True, re_calc = False):
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
			
	Note:
		the nearest neighbor is determined through the atomic configuration in
		initial_config_data
	"""
	
	# check if the strain_results_dict.pkl nn_results_dict.pkl file exists or not
	path_to_strain_results = path_to_test_dir + "/strain_results_dict.pkl"
	path_to_displacement = path_to_test_dir + "/displacement_results_dict.pkl"
	
	if re_calc is False:
		if os.path.exists(path_to_strain_results) and os.path.exists(path_to_displacement):
			print "atomic strain and displacement has already been calculated and saved in pkl file, skip"
			return (pickle.load(open(path_to_strain_results,'r')), pickle.load(open(path_to_displacement,'r')))
	
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
		init_atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"]==item])
		init_atom.box_dim = box_dim
		sad_atom = Atom.from_ds(saddle_config_data.loc[saddle_config_data["item"]==item])
		sad_atom.box_dim = box_dim
		disp_results[item] = Atom.distance_pbc(init_atom, sad_atom)
		
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
	Atom_ini_obj = Atom.from_ds(Atom_initial)
	
	
	NN_saddle = saddle_config_atom.loc[saddle_config_atom['item'] != atom_item]
	Atom_saddle = saddle_config_atom.loc[saddle_config_atom['item'] == atom_item]
	Atom_sad_obj = Atom.from_ds(Atom_saddle)
	
	Dim = 3
	V = np.zeros(shape=(Dim,Dim))
	W = np.zeros(shape=(Dim,Dim))
	
	for (index,atom_ini_NN) in NN_initial.iterrows():
		# d0_ji in pandas.Series
		d0_ji = Atom.from_ds(atom_ini_NN) - Atom_ini_obj
		d0_ji = Atom.to_list(d0_ji)
		atom_sad_NN = Atom.from_ds(NN_saddle.loc[NN_saddle["item"] == atom_ini_NN["item"]])
		# d_ji in pandas.Series
		d_ji = atom_sad_NN - Atom_sad_obj
		d_ji = Atom.to_list(d_ji)
		
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
