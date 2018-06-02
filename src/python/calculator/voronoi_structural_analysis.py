"""voronoi analysis module under periodic boudnary condition https://github.com/joe-jordan/pyvoro"""
import pyvoro
import os
import numpy as np
import pickle
import pandas as pd
import json
from collections import Counter
from util import operation_on_events

def run_all_tests_voronoi_calculator(path_to_data_dir,input_param):
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	
	box_range = input_param["box_range"]
	cut_off = input_param["cuf_off"]
	atom_list = input_param["atom_list"]
	periodic = input_param["periodic"]
	re_calc = input_param["re_calc"]
	
	operation = lambda x: single_event_voronoi_calculator(x, path_to_data_dir, box_range, cut_off, atom_list = atom_list, periodic = periodic, re_calc = re_calc)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc=num_of_proc)
	print "done voronoi index calculation for all interested tests!"

def single_event_voronoi_calculator(event_state, path_to_data_dir, box_range, cut_off, atom_list = None, periodic=[True,True,True], save_results = True, re_calc=False):
	"""
	this function calculates the voronoi index for all configurations in the test,
	on user specified atoms stored in atom_list
	"""
	# for local mode of voronoi calculation
	if type(atom_list) == dict:
	
	
	path_to_curr_result = path_to_test_dir + "/results"
	path_to_event_list = path_to_curr_result + "/selected_events.json"
	if not os.path.exists(path_to_curr_result):
		os.makedirs(path_to_curr_result)
	
	box_dim = [box_range[0][1] - box_range[0][0], box_range[1][1] - box_range[1][0], box_range[2][1] - box_range[2][0]]
	
	
	
	
	#escape when this test do not have log.file.1 file or configuration file
	try:
		event_list = event_selection(path_to_test_dir,box_dim,re_calc = re_calc)
	except IOError:
		return None
	# if no events has been selected in this test
	if event_list is None:
		return None
	
	# if final_selected_events.json exists, this means that event redudancy check
	# happens before calculating all strains
	# check if each event is in the final_selected_events.json file
	path_to_data_dir,test_id = os.path.split(path_to_test_dir)
	path_to_final_events = path_to_data_dir + "/final_selected_events.json"
	event_list = event_list.values()
	
	
	
	if save_results is True:
		print "begin saving voronoi results into json file"
		with open(path_to_voro_results, 'w') as f:
			json.dump(voronoi_index,f)
			f.close()


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
			
	
			
			
		
	
