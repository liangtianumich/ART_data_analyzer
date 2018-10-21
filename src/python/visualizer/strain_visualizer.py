"""
this strain visualization module plot the distribution of local atomic strains
in terms of atomic coordinates x,y,z in the simulation box
"""
import os
import json
import pickle
import pandas as pd
from util import Atom, data_dir_to_test_dir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from util import event_strain_disp,operation_on_events
from visualizer.general_visualizer import plot_histogram_3
# install mpl_toolkits module may need to add the __init__.py manually in
# site-packages/ to make it a package, if installed mpl_toolkits
# if mpl_tookits package not included along with matplotlib, sudo pip install -U matplotlib
# should do the trick, if not, i tried brew install geos, then sudo pip install - U matplotlib,
# it upgrades successfully, then sys.path updated and importlib.import_module('mpl_toolkits').__path__
# points to the symlink, brew remove python@2 --ignore-dependencies
# to remove the python2 dependencies created by brew to use system default python dependencies
def events_strain_visualization(path_to_data_dir, input_param):
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_strain_visualization(path_to_data_dir,x), num_of_proc)

def events_strain_visualization_old(path_to_data_dir, list_of_test_id):
	"""
	this function visualize/plot the strains for each events listed in tests in
	list_of_test_id
	
	plt.savefig
	if len(list_of_test_id) == 1:
		plt.show()
	"""
	test_id = ["test%s"%i for i in list_of_test_id]
	
	path_to_final_selected_events = path_to_data_dir + "final_selected_events.json"
	if os.path.exists(path_to_final_selected_events):
		final_selected_events = json.load(open(path_to_final_selected_events,"r"))
		final_interested_events = []
		for event in final_selected_events:
			if event[0] in test_id:
				final_interested_events.append(event)
		for event in final_interested_events:
			single_event_strain_visualization(path_to_data_dir, event)
	else:
		for test in list_of_test_id:
			path_to_test_dir = data_dir_to_test_dir(path_to_data_dir,test)
			path_to_test_result =  path_to_test_dir + "/results"
			path_to_event_list = path_to_test_result + "/selected_events.json"
			if os.path.exists(path_to_event_list):
				event_list = json.load(open(path_to_event_list,"r"))
				for value in event_list.values():
					event = ["test%s"%test,[value[0],value[1],value[2]]]
					single_event_strain_visualization(path_to_data_dir, event)
			else:
				print "skip current test:", "test%s"%test, "there is no selected events"	
	print "done plotting for the interested tests whose test_id is in the list",list_of_test_id

	
def single_event_strain_visualization(path_to_data_dir, event):
	"""
	this function plot the shear strain volumetric strain and displacement for a single event 
	"""	
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	path_to_test_result = path_to_test_dir + "/results"
	init, sad, fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_test_result + "/event_" + init + "_" + sad + "_" + fin
	path_to_init_sad = path_to_curr_event + "/init_sad"
	path_to_sad_fin = path_to_curr_event + "/sad_fin"
	path_to_init_fin = path_to_curr_event + "/init_fin"
	
	path_to_init_sad_strain = path_to_init_sad + "/strain_results_dict.pkl"
	path_to_init_sad_displacement = path_to_init_sad + "/displacement_results_dict.pkl"
	
	path_to_sad_fin_strain = path_to_sad_fin  + "/strain_results_dict.pkl"
	path_to_sad_fin_displacement = path_to_sad_fin + "/displacement_results_dict.pkl"
	
	path_to_init_fin_strain = path_to_init_fin + "/strain_results_dict.pkl"
	path_to_init_fin_displacement = path_to_init_fin + "/displacement_results_dict.pkl"
	
	init_sad_strain = pickle.load(open(path_to_init_sad_strain,'r'))
	init_sad_displacement = pickle.load(open(path_to_init_sad_displacement,'r'))
	
	sad_fin_strain = pickle.load(open(path_to_sad_fin_strain,'r'))
	sad_fin_displacement = pickle.load(open(path_to_sad_fin_displacement,'r'))
	
	init_fin_strain = pickle.load(open(path_to_init_fin_strain,'r'))
	init_fin_displacement = pickle.load(open(path_to_init_fin_displacement,'r'))
	
	init_sad_vol_strain, init_sad_shear_strain, init_sad_disp = event_strain_disp(init_sad_strain,init_sad_displacement)
	sad_fin_vol_strain, sad_fin_shear_strain, sad_fin_disp = event_strain_disp(sad_fin_strain, sad_fin_displacement)
	init_fin_vol_strain, init_fin_shear_strain, init_fin_disp = event_strain_disp(init_fin_strain, init_fin_displacement)
	
	path_to_init_sad_disp_strain = path_to_init_sad + '/disp_shear_strain.png'
	plot_2d_shear(path_to_init_sad_disp_strain,init_sad_disp,init_sad_shear_strain)
	
	path_to_init_sad_disp_vol_strain = path_to_init_sad + '/disp_vol_strain.png'
	plot_2d_vol(path_to_init_sad_disp_vol_strain,init_sad_disp,init_sad_vol_strain)
	
	
	path_to_sad_fin_disp_strain = path_to_sad_fin + '/disp_shear_strain.png'
	plot_2d_shear(path_to_sad_fin_disp_strain,sad_fin_disp,sad_fin_shear_strain)

	path_to_sad_fin_disp_vol_strain = path_to_sad_fin + '/disp_vol_strain.png'
	plot_2d_vol(path_to_sad_fin_disp_vol_strain,sad_fin_disp,sad_fin_vol_strain)
	
	path_to_init_fin_disp_strain = path_to_init_fin + '/disp_shear_strain.png'
	plot_2d_shear(path_to_init_fin_disp_strain,init_fin_disp,init_fin_shear_strain)
	
	path_to_init_fin_disp_vol_strain = path_to_init_fin + '/disp_vol_strain.png'
	plot_2d_vol(path_to_init_fin_disp_vol_strain,init_fin_disp,init_fin_vol_strain)
	
	
	plot_histogram_3(path_to_curr_event + "/disp_histogram.png", [init_sad_disp,sad_fin_disp,init_fin_disp])
	plot_histogram_3(path_to_curr_event + "/shear_strain_histogram.png", [init_sad_shear_strain,sad_fin_shear_strain,init_fin_shear_strain])
	plot_histogram_3(path_to_curr_event + "/vol_strain_histogram.png", [init_sad_vol_strain,sad_fin_vol_strain,init_fin_vol_strain])
	
	print "done plotting for the current event:"+ event[0] + "/event_" + init + "_" + sad + "_" + fin
			
def strain_events_stats_visualization(path_to_data_dir, input_param):
	"""
	use this new version of strain_events_stats_visualization need to reran the
	new strain_calc.py to overwrite events_stats.pkl with added event_state
	"""
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	all_events_results = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_strain_stats(path_to_data_dir,x), num_of_proc)
	
	
	disp_ave, disp_std, disp_max , disp_ave_2, disp_std_2, disp_max_2, disp_ave_3, disp_std_3, disp_max_3 = [], [], [], [], [], [], [], [], []
	
	shear_ave, shear_std, shear_max, shear_ave_2, shear_std_2, shear_max_2, shear_ave_3, shear_std_3, shear_max_3 = [], [], [], [], [], [], [], [], []
	
	vol_ave, vol_std, vol_max, vol_ave_2, vol_std_2, vol_max_2, vol_ave_3, vol_std_3, vol_max_3 = [], [], [], [], [], [], [], [], []
	
	for event_res in all_events_results:
		init_sad = event_res[0]
		sad_fin = event_res[1]
		init_fin = event_res[2]
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
	plot_histogram_3(path_to_data_dir+"/disp_ave.png", [disp_ave,disp_ave_2,disp_ave_3])
	plot_histogram_3(path_to_data_dir+"/disp_std.png", [disp_std,disp_std_2,disp_std_3])
	plot_histogram_3(path_to_data_dir+"/disp_max.png", [disp_max,disp_max_2,disp_max_3])
	
	plot_histogram_3(path_to_data_dir+"/shear_ave.png", [shear_ave,shear_ave_2,shear_ave_3])
	plot_histogram_3(path_to_data_dir+"/shear_std.png", [shear_std,shear_std_2,shear_std_3])
	plot_histogram_3(path_to_data_dir+"/shear_max.png", [shear_max,shear_max_2,shear_max_3])
	
	plot_histogram_3(path_to_data_dir+"/vol_ave.png", [vol_ave,vol_ave_2,vol_ave_3])
	plot_histogram_3(path_to_data_dir+"/vol_std.png", [vol_std,vol_std_2,vol_std_3])
	plot_histogram_3(path_to_data_dir+"/vol_max.png", [vol_max,vol_max_2,vol_max_3])	
	print "done plotting strain statistics for all interested tests!"
		
def single_event_strain_stats(path_to_data_dir,event):
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_curr_test = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	#path_to_curr_test = path_to_data_dir + event[0]
	path_to_curr_event = path_to_curr_test + "/results/events_stats.pkl"
	if os.path.exists(path_to_curr_event):
		result = pickle.load(open(path_to_curr_event,'r'))
		for event_res in result:
			# event_res[3] is the event state, if equal to curr event state: event
			if event_res[3] == event:
				init_sad = event_res[0]
				sad_fin = event_res[1]
				init_fin = event_res[2]
	return (init_sad, sad_fin, init_fin)	
	
def strain_events_stats_visualization_old(path_to_data_dir, list_of_test_id):
	"""
	this function visualize the strain statistics from tests listed in
	list_of_test_id, the user can customize the tests they want to know
	their statistics, among these tests, only the one has events_stat.pkl
	file will be used for plot
	
	path_to_data_dir: str
		path to the data dir
	list_of_tests: list
		a number list showing the test id number
		e.g. list [1,2] mean that visualizing test1, test2 data
		
	"""
	# set up the statistical quantities for the tests
	disp_ave, disp_std, disp_max , disp_ave_2, disp_std_2, disp_max_2, disp_ave_3, disp_std_3, disp_max_3 = [], [], [], [], [], [], [], [], []
	
	shear_ave, shear_std, shear_max, shear_ave_2, shear_std_2, shear_max_2, shear_ave_3, shear_std_3, shear_max_3 = [], [], [], [], [], [], [], [], []
	
	vol_ave, vol_std, vol_max, vol_ave_2, vol_std_2, vol_max_2, vol_ave_3, vol_std_3, vol_max_3 = [], [], [], [], [], [], [], [], []
	
	for i in list_of_test_id:
		path_to_curr_test = data_dir_to_test_dir(path_to_data_dir, i)
		path_to_curr_event = path_to_curr_test + "/results/events_stats.pkl"
		# skip the test who do not have events_stats.pkl in all tests specified in
		# list_of_test_id
		if os.path.exists(path_to_curr_event):
			result = pickle.load(open(path_to_curr_event,'r'))
			for event in result:
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
	
	plot_histogram_3(path_to_data_dir+"/disp_ave.png", [disp_ave,disp_ave_2,disp_ave_3])
	plot_histogram_3(path_to_data_dir+"/disp_std.png", [disp_std,disp_std_2,disp_std_3])
	plot_histogram_3(path_to_data_dir+"/disp_max.png", [disp_max,disp_max_2,disp_max_3])
	
	plot_histogram_3(path_to_data_dir+"/shear_ave.png", [shear_ave,shear_ave_2,shear_ave_3])
	plot_histogram_3(path_to_data_dir+"/shear_std.png", [shear_std,shear_std_2,shear_std_3])
	plot_histogram_3(path_to_data_dir+"/shear_max.png", [shear_max,shear_max_2,shear_max_3])
	
	plot_histogram_3(path_to_data_dir+"/vol_ave.png", [vol_ave,vol_ave_2,vol_ave_3])
	plot_histogram_3(path_to_data_dir+"/vol_std.png", [vol_std,vol_std_2,vol_std_3])
	plot_histogram_3(path_to_data_dir+"/vol_max.png", [vol_max,vol_max_2,vol_max_3])	
	print "done plotting strain statistics for all interested tests!"

def plot_2d_shear(path_to_dir, displacement, strain):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(displacement, strain,'ro',markersize=1.5)
	plt.xlabel('atomic displacement',fontsize=20)
	plt.ylabel('von Mises shear strain',fontsize=20)
	plt.savefig(path_to_dir)
	plt.close()

def plot_2d_vol(path_to_dir, displacement, strain):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(displacement, strain,'ro',markersize=1.5)
	plt.xlabel('atomic displacement',fontsize=20)
	plt.ylabel('volumetric strain',fontsize=20)
	plt.savefig(path_to_dir)
	plt.close()
	
def plot_strain_2d(initial_config_data, strain, normal_plane):
	"""
	
	this function plot the distributions of local atomic strains
	in 2d for von mises strain and hydrostatic invariants
	
	"""
	for key,atom_strain in strain.items():
		
		atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
		atom_strain[0]
		atom_strain[1]
		
