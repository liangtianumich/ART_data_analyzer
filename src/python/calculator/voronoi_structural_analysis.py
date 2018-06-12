"""voronoi analysis module under periodic boudnary condition https://github.com/joe-jordan/pyvoro"""
import pyvoro
import os
import numpy as np
import pickle
import pandas as pd
import json
from collections import Counter
from data_reader import *
from util import operation_on_events, event_local_atom_index, read_from_art_input_file
from visualizer.voronoi_visualizer import plot_voronoi_histogram_3, plot_dynamic_transition_matrix

# voronoi index classification from Evan Ma paper "Tuning order in disorder"
global ICO
ICO = [[0,4,4,0],[0,3,6,0],[0,2,8,0],[0,2,8,1],[0,0,12,0],[0,1,10,2],[0,0,12,2],[0,0,12,3],[0,0,12,4],[0,0,12,5]]
global ICO_LIKE
ICO_LIKE = [[0,5,2,1],[0,4,4,1],[0,3,6,1],[0,3,6,2],[0,2,8,2],[0,2,8,3],[0,1,10,3],[0,1,10,4],[0,1,10,5],[0,1,10,6],\
[0,6,0,2],[0,5,2,2],[0,4,4,2],[0,4,4,3],[0,3,6,3],[0,3,6,4],[0,2,8,4],[0,2,8,5],[0,2,8,6],[0,2,8,7]]
#global GUM
#GUM = [[0,6,0,3], [0,5,2,3], [0,5,2,4], [0,4,4,4], [0,4,4,5], [0,3,6,5], [0,3,6,6], [0,3,6,7], [0,3,6,8], \
#[0,6,0,4], [0,6,0,5], [0,5,2,5], [0,5,2,6], [0,4,4,6], [0,4,4,7], [0,4,4,8], [0,4,4,9]]

def run_all_tests_voronoi_calculator(path_to_data_dir, input_param):
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	
	box_range = input_param["box_range"]
	cut_off = input_param["cut_off"]
	atom_list = input_param["atom_list"]
	periodic = input_param["periodic"]
	re_calc = input_param["re_calc"]
	
	operation = lambda x: single_event_voronoi_calculator(x, path_to_data_dir, box_range, cut_off, atom_list = atom_list, periodic = periodic, re_calc = re_calc)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	print ("done voronoi index calculation for all interested tests!")

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
	
	for event_result in result_list:
		init_voronoi_class, sad_voronoi_class,fin_voronoi_class = event_result["init"], event_result["sad"], event_result["fin"]
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
		
	# begin calculate the probability for dynamic transition from init to sad
	init_total = total_init_ICO + total_init_ICO_LIKE + total_init_GUM
	sad_total = total_sad_ICO + total_sad_ICO_LIKE + total_sad_GUM
	init_ICO_pt = float(total_init_ICO)/init_total
	init_ICO_LIKE_pt = float(total_init_ICO_LIKE)/init_total
	init_GUM_pt = float(total_init_GUM)/init_total
	
	p11_0 = 1.0/3 
	p12_0 = p11_0
	p13_0 = p11_0
	p21_0 = 1.0/3
	p22_0 = p21_0
	p23_0 = p21_0
	p31_0 = 1.0/3
	p32_0 = p31_0
	p33_0 = p31_0
	
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
	print p
	print c_matrix
	
	path_to_image = path_to_data_dir + "/dynamic_transition_matrix_all_events.png"
	plot_dynamic_transition_matrix(path_to_image, c_matrix)
	
	print ("done voronoi index classification for all interested tests!")
	
def single_event_voronoi_classifier(event_state, path_to_data_dir):
	"""
	this function load the calculated voronoi index results file voronoi_index_results.json
	classify the vornoi index into ICO, ICO_LIKE, GUM according to the criterion
	defined at the top level of the module
	"""
	path_to_test_dir = path_to_data_dir + event_state[0]
	path_to_curr_result = path_to_test_dir + "/results"
	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	event_str = event_state[0] + "/event_" + init + "_" + sad + "_" + fin
	if not os.path.exists(path_to_curr_event):
		raise Exception("the voronoi index has not been calculated for the event %s"%event_str)
	
	path_to_voro_results = path_to_curr_event + "/voronoi_index_results.json"
	if not os.path.exists(path_to_voro_results):
		raise Exception("the voronoi index has not been calculated for the event %s"%event_str)
	
	voronoi_index = json.load(open(path_to_voro_results,"r"))
	# classify voronoi index
	init_voronoi_class = classify_voronoi_index(voronoi_index["init"])
	sad_voronoi_class = classify_voronoi_index(voronoi_index["sad"])
	fin_voronoi_class = classify_voronoi_index(voronoi_index["fin"])
	
	voronoi_class = {"init": init_voronoi_class, "sad": sad_voronoi_class, "fin": fin_voronoi_class}
	
	# do visualization of the fraction of voronoi class
	path_to_image = path_to_curr_event + "/voronoi_hist.png"
	plot_voronoi_histogram_3(path_to_image, [init_voronoi_class,sad_voronoi_class,fin_voronoi_class])
	
	return voronoi_class
	
def single_event_voronoi_calculator(event_state, path_to_data_dir, box_range, cut_off, atom_list = None,max_edge_count=8, periodic = [True,True,True], save_results = True, re_calc = False):
	"""
	this function calculates the voronoi index of user specified atoms stored in atom_list
	for all configurations (init,sad,fin) in a single event that are specified in the event state
	Input:
		event_state: a list
			a list with the 1st element being the test_id, e.g. test1
			the 2nd element being a list containing the string of init, sad, fin
			configuration file str, e.g. [min1000,sad1001,min1001]
	"""
	
			
	path_to_test_dir = path_to_data_dir + event_state[0]
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
	if re_calc is False:
		if os.path.exists(path_to_voro_results):
			return json.load(open(path_to_voro_results,"r"))
	
	path_to_file_ini = path_to_test_dir + '/' + init + ".dump"
	path_to_file_sad = path_to_test_dir + '/' + sad + ".dump"
	path_to_file_fin = path_to_test_dir + '/' + fin + ".dump"
		
	initial_config_data = read_data_from_file(path_to_file_ini)
	saddle_config_data = read_data_from_file(path_to_file_sad)
	final_config_data = read_data_from_file(path_to_file_fin)
	
	box_dim = [box_range[0][1] - box_range[0][0], box_range[1][1] - box_range[1][0], box_range[2][1] - box_range[2][0]]
	
	path_to_local_atom_index = path_to_curr_event + "/local_atoms_index.json"
	
	if atom_list is None:
		if os.path.exists(path_to_local_atom_index):
			# for local mode of voronoi calculation
			print ("\n starting local mode voronoi calculations")
			local_atom_list = json.load(open(path_to_local_atom_index,"r"))
			atom_list = [atom + 1 for atom in local_atom_list]
		else:
			atom_list = (initial_config_data["item"]).tolist()

	init_voronoi_index = single_config_voronoi_calculator(initial_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
	sad_voronoi_index = single_config_voronoi_calculator(saddle_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
	fin_voronoi_index = single_config_voronoi_calculator(final_config_data, box_range, cut_off, atom_list=atom_list, max_edge_count = max_edge_count, periodic=periodic)
	
	voronoi_index = {"init":init_voronoi_index, "sad":sad_voronoi_index, "fin":fin_voronoi_index}
	
	if save_results is True:
		print ("begin saving voronoi results into json file")
		with open(path_to_voro_results, 'w') as f:
			json.dump(voronoi_index,f)
			f.close()
	return voronoi_index


def single_config_voronoi_calculator(config, box_range, cut_off, atom_list=None, max_edge_count=8, periodic=[True,True,True]):
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
	
	return:
		voronoi_index: list of lists
			voronoi_index for each atom in atom_list
	"""
	
	# read the configuration data
	points = (config[['x','y','z']].values).tolist()
	
	if atom_list is None:
		atom_list = (config["item"].values).tolist()
	int_points_index = [atom-1 for atom in atom_list]
	
	dispersion = cut_off
	
	all_results = pyvoro.compute_voronoi(points, box_range, dispersion, periodic=periodic)
	
	results = [all_results[i] for i in int_points_index]
	
	voronoi_index = count_faces(results, max_edge_count)
	
	return voronoi_index

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
		
def count_faces(results, max_edge_count=8):
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
	return all_voronoi_index
			
	
			
			
		
	
