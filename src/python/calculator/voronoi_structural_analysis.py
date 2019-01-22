"""voronoi analysis module under periodic boudnary condition https://github.com/joe-jordan/pyvoro"""
import tess
import pyvoro
import os
import numpy as np
import pickle
import pandas as pd
import json
from collections import Counter
from data_reader import *
from util import operation_on_events, event_local_atom_index, read_from_art_input_file, data_dir_to_test_dir
from visualizer.voronoi_visualizer import plot_voronoi_histogram_3, plot_dynamic_transition_matrix
# voronoi index classification from Evan Ma paper "Tuning order in disorder"
global ICO
ICO = [[0,6,0,0],[0,5,2,0],[0,4,4,0],[0,3,6,0],[0,2,8,0],[0,2,8,1],[0,0,12,0],[0,1,10,2],[0,0,12,2],[0,0,12,3],[0,0,12,4],[0,0,12,5]]
global ICO_LIKE
ICO_LIKE = [[0,6,0,1],[0,5,2,1],[0,4,4,1],[0,3,6,1],[0,3,6,2],[0,2,8,2],[0,2,8,3],[0,1,10,3],[0,1,10,4],[0,1,10,5],[0,1,10,6],\
[0,6,0,2],[0,5,2,2],[0,4,4,2],[0,4,4,3],[0,3,6,3],[0,3,6,4],[0,2,8,4],[0,2,8,5],[0,2,8,6],[0,2,8,7]]
#global GUM
#GUM = [[0,6,0,3], [0,5,2,3], [0,5,2,4], [0,4,4,4], [0,4,4,5], [0,3,6,5], [0,3,6,6], [0,3,6,7], [0,3,6,8], \
#[0,6,0,4], [0,6,0,5], [0,5,2,5], [0,5,2,6], [0,4,4,6], [0,4,4,7], [0,4,4,8], [0,4,4,9]]

def run_all_tests_voronoi_calculator(path_to_data_dir, input_param, return_volume = False):
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	
	box_range = input_param["box_range"]
	cut_off = input_param["cut_off"]
	atom_list = input_param["atom_list"]
	periodic = input_param["periodic"]
	re_calc = input_param["re_calc"]
	
	operation = lambda x: single_event_voronoi_calculator(x, path_to_data_dir, box_range, cut_off, atom_list = atom_list, periodic = periodic, re_calc = re_calc, return_volume = return_volume)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	print ("done voronoi cell calculations for all interested tests!")

def run_all_tests_voronoi_classifier(path_to_data_dir, input_param):
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	
	operation = lambda x: single_event_voronoi_classifier(x, path_to_data_dir)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	total_init_ICO, total_init_ICO_LIKE, total_init_GUM = 0,0,0
	total_sad_ICO, total_sad_ICO_LIKE, total_sad_GUM = 0,0,0
	total_fin_ICO, total_fin_ICO_LIKE, total_fin_GUM = 0,0,0
	ICO_to_ICO = 0
	ICO_to_ICO_LIKE = 0
	ICO_to_GUM = 0
	ICO_LIKE_to_ICO = 0
	ICO_LIKE_to_ICO_LIKE = 0
	ICO_LIKE_to_GUM = 0
	GUM_to_ICO = 0
	GUM_to_ICO_LIKE = 0
	GUM_to_GUM = 0
	init_voronoi_class_tol, sad_voronoi_class_tol, fin_voronoi_class_tol = [], [], []
	for event_result in result_list:
		if event_result is None:
			continue
		init_voronoi_class, sad_voronoi_class,fin_voronoi_class = event_result["init"], event_result["sad"], event_result["fin"]
		
		init_voronoi_class_tol.extend(init_voronoi_class)
		sad_voronoi_class_tol.extend(sad_voronoi_class)
		fin_voronoi_class_tol.extend(fin_voronoi_class)
		
		atom_index = range(len(init_voronoi_class))
		for atom_id in atom_index:
			if init_voronoi_class[atom_id] == 0:
				if sad_voronoi_class[atom_id] == 0:
					ICO_to_ICO = ICO_to_ICO + 1
				elif sad_voronoi_class[atom_id] == 1:
					ICO_to_ICO_LIKE = ICO_to_ICO_LIKE + 1
				elif sad_voronoi_class[atom_id] == 2:
					ICO_to_GUM = ICO_to_GUM + 1

			if init_voronoi_class[atom_id] == 1:
				if sad_voronoi_class[atom_id] == 0:
					ICO_LIKE_to_ICO = ICO_LIKE_to_ICO + 1
				elif sad_voronoi_class[atom_id] == 1:
					ICO_LIKE_to_ICO_LIKE = ICO_LIKE_to_ICO_LIKE + 1
				elif sad_voronoi_class[atom_id] == 2:
					ICO_LIKE_to_GUM = ICO_LIKE_to_GUM + 1
			
			if init_voronoi_class[atom_id] == 2:
				if sad_voronoi_class[atom_id] == 0:
					GUM_to_ICO = GUM_to_ICO + 1
				elif sad_voronoi_class[atom_id] == 1:
					GUM_to_ICO_LIKE = GUM_to_ICO_LIKE + 1
				elif sad_voronoi_class[atom_id] == 2:
					GUM_to_GUM = GUM_to_GUM + 1
			
		init_count = Counter(init_voronoi_class)
		sad_count = Counter(sad_voronoi_class)
		fin_count = Counter(fin_voronoi_class)
		
		total_init_ICO = total_init_ICO + init_count[0]
		total_init_ICO_LIKE = total_init_ICO_LIKE + init_count[1]
		total_init_GUM = total_init_GUM + init_count[2]
		
		total_sad_ICO = total_sad_ICO + sad_count[0]
		total_sad_ICO_LIKE = total_sad_ICO_LIKE + sad_count[1]
		total_sad_GUM = total_sad_GUM + sad_count[2]
		
		total_fin_ICO = total_fin_ICO + fin_count[0]
		total_fin_ICO_LIKE = total_fin_ICO_LIKE + fin_count[1]
		total_fin_GUM = total_fin_GUM + fin_count[2]
		# work on more statistics if necessary
		
	
	init_total = total_init_ICO + total_init_ICO_LIKE + total_init_GUM
	sad_total = total_sad_ICO + total_sad_ICO_LIKE + total_sad_GUM
	fin_total = total_fin_ICO + total_fin_ICO_LIKE + total_fin_GUM
	
	init_ICO_pt = float(total_init_ICO)/init_total
	init_ICO_LIKE_pt = float(total_init_ICO_LIKE)/init_total
	init_GUM_pt = float(total_init_GUM)/init_total
	init_pt = [init_ICO_pt, init_ICO_LIKE_pt, init_GUM_pt]
	
	sad_ICO_pt = float(total_sad_ICO)/sad_total
	sad_ICO_LIKE_pt = float(total_sad_ICO_LIKE)/sad_total
	sad_GUM_pt = float(total_sad_GUM)/sad_total
	sad_pt = [sad_ICO_pt, sad_ICO_LIKE_pt, sad_GUM_pt]
	
	fin_ICO_pt = float(total_fin_ICO)/fin_total
	fin_ICO_LIKE_pt = float(total_fin_ICO_LIKE)/fin_total
	fin_GUM_pt = float(total_fin_GUM)/fin_total
	fin_pt = [fin_ICO_pt, fin_ICO_LIKE_pt, fin_GUM_pt]
	path_to_voro_class_pt = path_to_data_dir + "/voronoi_class_fraction_all_events.png"
	print "All filtered events in list_of_test_id:"
	print "initial state ICO, ICO-like, GUM fraction is:", init_pt
	print "sadlle state ICO, ICO-like, GUM fraction is:", sad_pt
	print "final state ICO, ICO-like, GUM fraction is:", fin_pt
	plot_voronoi_histogram_3(path_to_voro_class_pt, [init_voronoi_class_tol, sad_voronoi_class_tol, fin_voronoi_class_tol])
	
	# begin calculate the probability for dynamic transition from init to sad
	p11_0 = 1.0/3 
	p12_0 = p11_0
	p13_0 = p11_0
	p21_0 = 1.0/3
	p22_0 = p21_0
	p23_0 = p21_0
	p31_0 = 1.0/3
	p32_0 = p31_0
	p33_0 = p31_0
	
	if total_init_ICO == 0 or total_init_ICO_LIKE == 0 or total_init_GUM == 0:
		print "total number of ICO is:", total_init_ICO
		print "total number of ICO-LIKE is:", total_init_ICO_LIKE
		print "total number of GUM is:", total_init_GUM
		print "Can not calculate the dynamic transition probability matrix since either total number of ICO or ICO-LIKE or GUM is zero!"
		return
	
	p11 = float(ICO_to_ICO)/total_init_ICO
	p12 = float(ICO_to_ICO_LIKE)/total_init_ICO
	p13 = float(ICO_to_GUM)/total_init_ICO
	p21 = float(ICO_LIKE_to_ICO)/total_init_ICO_LIKE
	p22 = float(ICO_LIKE_to_ICO_LIKE)/total_init_ICO_LIKE
	p23 = float(ICO_LIKE_to_GUM)/total_init_ICO_LIKE
	p31 = float(GUM_to_ICO)/total_init_GUM
	p32 = float(GUM_to_ICO_LIKE)/total_init_GUM
	p33 = float(GUM_to_GUM)/total_init_GUM
	
	p = np.array([[p11,p12,p13],[p21,p22,p23],[p31,p32,p33]])
	p_0 = np.array([[p11_0,p12_0,p13_0],[p21_0,p22_0,p23_0],[p31_0,p32_0,p33_0]])
	
	c_matrix = p/p_0 - 1
	print "Probability matrix:", p
	print "Normalized probability matrix:", c_matrix
	
	path_to_image = path_to_data_dir + "/dynamic_transition_probability_matrix_all_events.png"
	plot_dynamic_transition_matrix(path_to_image, p)
	
	print ("done voronoi index classification for all interested tests!")
	
def single_event_voronoi_classifier(event_state, path_to_data_dir):
	"""
	this function load the calculated voronoi index results file voronoi_index_results.json
	classify the vornoi index into ICO, ICO_LIKE, GUM according to the criterion
	defined at the top level of the module
	"""
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	#path_to_test_dir = path_to_data_dir + event_state[0]
	path_to_curr_result = path_to_test_dir + "/results"
	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	event_str = event_state[0] + "/event_" + init + "_" + sad + "_" + fin
	
	path_to_voro_results = path_to_curr_event + "/voronoi_index_results.json"
	if not os.path.exists(path_to_voro_results):
		print("the voronoi index has not been calculated for the event %s"%event_str)
		return None
	
	voronoi_index = json.load(open(path_to_voro_results,"r"))
	if voronoi_index["init"] == []:
		return None
	# classify voronoi index
	init_voronoi_class = classify_voronoi_index(voronoi_index["init"])
	sad_voronoi_class = classify_voronoi_index(voronoi_index["sad"])
	fin_voronoi_class = classify_voronoi_index(voronoi_index["fin"])
	
	voronoi_class = {"init": init_voronoi_class, "sad": sad_voronoi_class, "fin": fin_voronoi_class}
	
	# do visualization of the fraction of voronoi class
	path_to_image = path_to_curr_event + "/voronoi_hist.png"
	plot_voronoi_histogram_3(path_to_image, [init_voronoi_class,sad_voronoi_class,fin_voronoi_class])
	
	return voronoi_class
	
def single_event_voronoi_calculator(event_state, path_to_data_dir, box_range, cut_off, atom_list = None,max_edge_count=8, periodic = [True,True,True], save_results = True, re_calc = False, return_volume=False):
	"""
	this function calculates the voronoi index of user specified atoms stored in atom_list
	for all configurations (init,sad,fin) in a single event that are specified in the event state
	Input:
		event_state: a list
			a list with the 1st element being the test_id, e.g. test1
			the 2nd element being a list containing the string of init, sad, fin
			configuration file str, e.g. [min1000,sad1001,min1001]
	"""
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
			
	#path_to_test_dir = path_to_data_dir + event_state[0]
	path_to_curr_result = path_to_test_dir + "/results"
	# redundant, if exists, may cause race condition when os.makedirs act on same dir (leaf)
	# for parallel processing, in python3, it can be avoided by adding exist_ok=True
	# in python2 here, rewrite code to ensure os.makedirs act on different dir (leaf)
	# or using try except in 
	# https://stackoverflow.com/questions/12468022/python-fileexists-error-when-making-directory
	# https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
	#if not os.path.exists(path_to_curr_result):
	#	os.makedirs(path_to_curr_result)
	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	if not os.path.exists(path_to_curr_event):
		os.makedirs(path_to_curr_event)
	
	path_to_voro_results = path_to_curr_event + "/voronoi_index_results.json"
	path_to_volume_results = path_to_curr_event + "/voronoi_volume_results.json"
	
	if re_calc is False:
		if return_volume is True:
			if os.path.exists(path_to_voro_results) and os.path.exists(path_to_volume_results):
				return (json.load(open(path_to_voro_results,"r")), json.load(open(path_to_volume_results,"r")))
			else:
				print "return volume is True, no existing calculations, begin new calculations!"
				#raise Exception("User specify return_volume of voronoi cell is True, however, either voronoi_index_results.json or voronoi_volume_results.json does not exists!")
		elif os.path.exists(path_to_voro_results):
			return json.load(open(path_to_voro_results,"r"))
		else:
			print "begin calculating voronoi index"
	else:
		print "re_calculating"
	path_to_file_ini = path_to_test_dir + '/' + init + ".dump"
	path_to_file_sad = path_to_test_dir + '/' + sad + ".dump"
	path_to_file_fin = path_to_test_dir + '/' + fin + ".dump"

	initial_config_data = read_data_from_file(path_to_file_ini)
	saddle_config_data = read_data_from_file(path_to_file_sad)
	final_config_data = read_data_from_file(path_to_file_fin)
	
	box_dim = [box_range[0][1] - box_range[0][0], box_range[1][1] - box_range[1][0], box_range[2][1] - box_range[2][0]]
	
	path_to_local_atom_index = path_to_curr_event + "/local_atoms_index.json"
	path_to_initial_atom_index = path_to_curr_event + "/initial_cluster_atoms_index.json"
	path_to_central_atom_index = path_to_curr_event + "/central_atom_index.json"
	path_to_max_disp_atom_index = path_to_curr_event + "/max_disp_atom_index.json"
	
	if atom_list is None:
		atom_list = (initial_config_data["item"]).tolist()
	
	if atom_list == "local":
		if os.path.exists(path_to_local_atom_index):
			# for local mode of voronoi calculation
			print ("\n starting local mode voronoi calculations for local atoms during init to sad process")
			local_atom_list = json.load(open(path_to_local_atom_index,"r"))
			if local_atom_list["init_sad"] == []:
				print "No local atoms during init to sad, return None"
				return None
			#atom_list = [atom + 1 for atom in local_atom_list["init_sad"]]
			atom_list = local_atom_list["init_sad"]
		else:
			raise Exception("local_atoms_index.json does not exist in %s"%path_to_curr_event)
	elif atom_list == "initial":
		if os.path.exists(path_to_initial_atom_index):
			print ("\n starting initial triggered cluster atoms voronoi calculations")
			initial_atom_list = json.load(open(path_to_initial_atom_index,"r"))
			if initial_atom_list == []:
				return None
			atom_list = initial_atom_list
		else:
			raise Exception("initial_cluster_atoms_index.json does not exist in %s"%path_to_curr_event)
	elif atom_list == "central":
		if os.path.exists(path_to_central_atom_index):
			print ("\n starting triggered central atom voronoi calculations")
			central_atom_list = json.load(open(path_to_central_atom_index,"r"))
			if central_atom_list == []:
				return None
			atom_list = central_atom_list
		else:
			raise Exception("central_atom_index.json does not exist in %s"%path_to_curr_event)
	elif atom_list == "max_disp":
		if os.path.exists(path_to_max_disp_atom_index):
			print ("\n starting the max disp atom voronoi calculations")
			max_disp_atom_list = json.load(open(path_to_max_disp_atom_index,"r"))
			if max_disp_atom_list == []:
				return None
			atom_list = max_disp_atom_list
		else:
			raise Exception("max_disp_atom_index.json does not exist in %s"%path_to_curr_event)

	if return_volume is True:
		init_voronoi_index, init_volumes = single_config_voronoi_calculator(initial_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic, return_volume=True)
		sad_voronoi_index, sad_volumes = single_config_voronoi_calculator(saddle_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic, return_volume=True)
		fin_voronoi_index, fin_volumes = single_config_voronoi_calculator(final_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic,return_volume=True)
		voronoi_volume = {"init":init_volumes, "sad":sad_volumes, "fin":fin_volumes}
	else:
		init_voronoi_index = single_config_voronoi_calculator(initial_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
		sad_voronoi_index = single_config_voronoi_calculator(saddle_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
		fin_voronoi_index = single_config_voronoi_calculator(final_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
	
	voronoi_index = {"init":init_voronoi_index, "sad":sad_voronoi_index, "fin":fin_voronoi_index}
	
	if save_results is True:
		print ("begin saving voronoi results into json file")
		with open(path_to_voro_results, 'w') as f:
			json.dump(voronoi_index,f)
			f.close()
		if return_volume is True:
			with open(path_to_volume_results, 'w') as f:
				json.dump(voronoi_volume,f)
				f.close()
	
	if return_volume is True:
		return (voronoi_index, voronoi_volume)
	else:
		return voronoi_index

def single_config_voronoi_calculator(config, box_range, cut_off, atom_list=None, max_edge_count=8, periodic=[True,True,True],return_volume=False, tool="pyvoro"):
	"""
	this function calculates the voronoi index for atoms in atom_list
	
	config: pandas.Dataframe
		configuration data stored in pandas.Dataframe
	
	atom_list: list or dict
		list of atom item id or dict for local mode
	
	box_range: list
		list of simulation box range in x, y, z direction
	
	cut_off: float
		max distance between two points that might be adjacent, which sets
		voro++ block sizes
	
	max_edge_count: int
		an integer that specify the maximum number of edges that will be trucated
		in voronoi index, this can be specified from user's experience for their problems
		or calculated from the maximum of all max_edge_count for all particles
	
	return_volume: Boolean
		if True, return the voronoi cell volume of each atomic particle along with their voronoi indexes

	return:
		voronoi_index: list of lists
			voronoi_index for each atom in atom_list
	"""
	
	# read the configuration data
	# points = (config[['x','y','z']].values).tolist()
	
	if atom_list is None:
		atom_list = (config["item"].values).tolist()
	
	int_points = config.loc[config['item'].isin(atom_list)]
	
	# order int_points based on the appearance of atom in atom_list
	df = pd.DataFrame()
	for atom in atom_list:
		df = df.append(int_points.loc[int_points['item'] == atom])
	int_points = df 
		
	#int_points_item_id = (int_points["item"].values).tolist()
	#int_points_index = [atom-1 for atom in atom_list]
	
	dispersion = cut_off
	
	box_dim = [box_range[0][1] - box_range[0][0], box_range[1][1] - box_range[1][0], box_range[2][1] - box_range[2][0]]
	
	int_voro_results = []
	for index,point in int_points.iterrows():
		# calculate NN of this point
		# result_group = result_tree.query_ball_point(point[['x','y','z']].values, cut_off * 2.0)
		[x_c, y_c, z_c] = point[['x','y','z']].values
				
		config_corr = config[config.item != point['item']]
		
		NN_df = NN_cube_pbc(config_corr,[x_c, y_c, z_c], box_range, cut_off)
		
		#NN_df = config.iloc[result_group]
		all_df = NN_df.append(point,ignore_index=True)
		all_points = all_df[['x','y','z']].values		
		# let this NN with this point to calculate voronoi indexes results
		if tool == "tess":
			limits = (box_range[0][0],box_range[1][0],box_range[2][0]),(box_range[0][1],box_range[1][1], box_range[2][1])
			results = tess.Container(all_points, limits, periodic=tuple(periodic))
		else:
			results = pyvoro.compute_voronoi(all_points, box_range, dispersion, periodic=periodic)
		
		# append the voronoi index of only this point
		curr_voro = results[-1]
		int_voro_results.append(curr_voro)

	if return_volume is True:
		voronoi_index, volumes = count_faces(int_voro_results, max_edge_count, True, tool=tool)
		return (voronoi_index, volumes)
	else:
		voronoi_index = count_faces(int_voro_results, max_edge_count, False,tool=tool)
		return voronoi_index

def NN_cube_pbc(config,point, box_range, cut_off):
	"""
	this function finds the nearest neighbor points in terms of cube around point
	within a cut_off in each dimension x,y,z direction
	
	Input:
		config: pandas.DataFrame
			all data points to search for NN
		point: a list
			a list contains coordinates in x,y,z direction
		box_range: a list of 3 lists
			e.g. [[x_low, x_high],[y_low, y_high], [z_low, z_high]]
		cuf_off: float
			the cut-off perpendicular distance from center to cube edge
	Return:
		config_z: pandas.DataFrame
			all NN data points within the cube considering periodic boundary
			condition
	"""
	[x,y,z] = point
	x_range = [x - cut_off * 1.1, x + cut_off * 1.1]
	y_range = [y - cut_off * 1.1, y + cut_off * 1.1]
	z_range = [z - cut_off * 1.1, z + cut_off * 1.1]
	[box_range_x, box_range_y, box_range_z] = [box_range[0],box_range[1],box_range[2]]
	config_x = df_range_pbc('x',x_range, box_range_x, config)
	config_y = df_range_pbc('y',y_range, box_range_y, config_x)
	config_z = df_range_pbc('z',z_range, box_range_z, config_y)
	return config_z

def df_range_pbc(axis,raw_range, dim_range, config):
	"""
	this function select a subset of data points in config based on the
	required range specified in raw_range in specified axis considering
	the periodic boundary condition.
	Input: 
		axis: str,
			e.g. 'x'
		
		raw_range: list
			e.g. [1.234, 3.769]
		
		dim_range: list
			e.g. [0.234, 35.978]
		
		config: pandas.DataFrame
			all data points
	
	Return:
		a pandas.DataFrame whose data points are within the range along specified axis
		considering pbc
	"""
	if raw_range[0] >= dim_range[0] and raw_range[1] <= dim_range[1]:
		# raw_range
		return config.loc[(config[axis] >= raw_range[0]) & (config[axis] <= raw_range[1])]
	elif raw_range[0] <= dim_range[0] and raw_range[1] <= dim_range[1]:
		#([dim_range[0], raw_range[1]], [raw_range[0] + dim_range[1] - dim_range[0],dim_range[1]])
		return config.loc[((config[axis] >= dim_range[0]) & (config[axis] <= raw_range[1])) | ((config[axis] >= raw_range[0] + dim_range[1] - dim_range[0]) & (config[axis] <= dim_range[1]))]
	elif raw_range[0] <= dim_range[1] and raw_range[1] >= dim_range[1]:
		#return ([raw_range[0], dim_range[1]], [dim_range[0], raw_range[1] - dim_range[1] + dim_range[0]])
		return config.loc[((config[axis] >= raw_range[0]) & (config[axis] <= dim_range[1])) | ((config[axis] >= dim_range[0]) & (config[axis] <= raw_range[1] - dim_range[1] + dim_range[0]))]
	else:
		raise Exception("either specified data range or simulation box range is wrong!")

def classify_voronoi_index(list_of_voronoi_index):
	"""
	this function classify a list of voronoi index into a list of
	ico, ico-like, GUMs
	"""
	list_of_voronoi_class = []
	for x in list_of_voronoi_index:
		if len(x) < 6:
			raise Exception("can not classify voronoi index vector whose length is less than 6")
		else:
			truncated_x = x[2:6]
		
		if truncated_x in ICO:
			#list_of_voronoi_class.append('ico')
			list_of_voronoi_class.append(0)
		elif truncated_x in ICO_LIKE:
			#list_of_voronoi_class.append('ico_like')
			list_of_voronoi_class.append(1)
		else:
			#list_of_voronoi_class.append('GUM')
			list_of_voronoi_class.append(2)
	return list_of_voronoi_class			

def count_faces(results, max_edge_count=8, return_volume=False, tool="tess"):
	if tool == "tess":
		return tess_count_faces(results, max_edge_count=max_edge_count, return_volume=return_volume)
	else:
		return pyvoro_count_faces(results, max_edge_count=max_edge_count, return_volume=return_volume)
		
def tess_count_faces(results, max_edge_count=8, return_volume=False):
	"""
	input:
		results: a list
			a list contains the voronoi cell information for each particle in the box
		max_edge_count: int
			an integer that specify the maximum number of edges that will be trucated
			in voronoi index, this can be specified from user's experience for their problems
			or calculated from the maximum of all max_edge_count for all particles
	return:
		all_voronoi_index: list
			a list contains the voronoi index vector for each particle.
			the voronoi index vector will be truncated at max_edge_count
	Note:
		the voronoi index vector here have a full representation, with
		the first two elements will always be zero.
	"""
	all_voronoi_index = []
	
	if return_volume is True:
		volumes = []
		for point in results:
			volumes.append(point.volume())

	for point in results:
		face_vertices = point.face_vertices()
		voronoi_index = []
		edges = []
		for face in face_vertices:
			edges.append(len(face))
		edges_count = Counter(edges)
		
		for i in xrange(max_edge_count):
			if i+1 in edges_count:
				voronoi_index.append(edges_count[i+1])
			else:
				voronoi_index.append(0)
		all_voronoi_index.append(voronoi_index)
	
	if return_volume is True:
		return (all_voronoi_index, volumes)
	else:
		return all_voronoi_index

def pyvoro_count_faces(results, max_edge_count=8, return_volume=False):
	"""
	input:
		results: a list
			a list contains the voronoi cell information for each particle in the box
		max_edge_count: int
			an integer that specify the maximum number of edges that will be trucated
			in voronoi index, this can be specified from user's experience for their problems
			or calculated from the maximum of all max_edge_count for all particles
	return:
		all_voronoi_index: list
			a list contains the voronoi index vector for each particle.
			the voronoi index vector will be truncated at max_edge_count
	Note:
		the voronoi index vector here have a full representation, with
		the first two elements will always be zero.
	"""
	all_voronoi_index = []
	
	if return_volume is True:
		volumes = []
		for point in results:
			volumes.append(point["volume"])
		
	for point in results:
		voronoi_index = []
		edges = []
		for face in point["faces"]:
			edges.append(len(face['vertices']))
		edges_count = Counter(edges)
		
		for i in xrange(max_edge_count):
			if i+1 in edges_count:
				voronoi_index.append(edges_count[i+1])
			else:
				voronoi_index.append(0)
		all_voronoi_index.append(voronoi_index)
	
	if return_volume is True:
		return (all_voronoi_index, volumes)
	else:
		return all_voronoi_index
			
	
			
			
		
	
