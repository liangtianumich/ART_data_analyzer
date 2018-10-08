"""
this module implements correlation model
"""
import os
import pickle
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
from sklearn import linear_model
from sklearn.svm import LinearSVC, LinearSVR
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import NearestNeighbors
from visualizer.general_visualizer import plot_histogram_3, plot_2d
from util import operation_on_events, Configuration, state_energy_barrier, data_dir_to_test_dir, read_from_art_log_file, read_cluster_radius_from_bart
from data_reader import read_data_from_file
from periodic_kdtree import PeriodicCKDTree

def shear_strain_vol_strain_cluster_all_events(path_to_data_dir, input_param, save_results=True):
	"""
	this function find cluster averaged shear strain vs cluster averaged volumetric strain 
	for all final selected events in list_of_test_id
	"""
	path_to_shear_vol_strain_init_sad_png = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_init_sad_all_events.png")
	#path_to_shear_vol_strain_init_fin_png = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_init_fin_all_events.png")
	#path_to_shear_vol_strain_sad_fin_png = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_sad_fin_all_events.png")
	
	path_to_shear_vol_init_sad = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_init_sad_all_events.json")
	#path_to_shear_vol_sad_fin = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_sad_fin_all_events.json")
	#path_to_shear_vol_init_fin = os.path.join(path_to_data_dir, "shear_vol_strain_cluster_init_fin_all_events.json")
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_cluster_averaged_shear_vol_strain(x, path_to_data_dir),num_of_proc)
	
	event_shear_strain, event_vol_strain, event_disp = [], [], []
	for result in result_list:
		event_shear_strain.append(result[0])
		event_vol_strain.append(result[1])
		#event_disp.append(result[2])
	
	plot_2d(path_to_shear_vol_strain_init_sad_png, event_shear_strain, event_vol_strain, "local event averaged shear strain","local event averaged volumetric strain")
	
	if save_results is True:
		with open(path_to_shear_vol_init_sad, 'w') as f:
			json.dump(result_list, f)
			f.close()
	print "done finding cluster averaged shear strain vs cluster averaged volumetric strain for all final selected events in list_of_test_id!"
	
def	single_event_cluster_averaged_shear_vol_strain(event, path_to_data_dir):
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	path_to_curr_result = path_to_test_dir + "/results"
	init,sad,fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	print "path_to_current_event:", path_to_curr_event
	path_to_local_atom_index = path_to_curr_event + "/local_atoms_index.json"
	if os.path.exists(path_to_local_atom_index):
		local_atoms_index = json.load(open(path_to_local_atom_index),'r')
	else:
		raise Exception("indexes of cluster local atoms has not been determined, please find the cluster local atoms first")
	
	path_to_init_sad = path_to_curr_event + "/init_sad"
	path_to_all_strain_results_init_sad = path_to_init_sad + "/strain_results_dict.pkl"
	path_to_all_displacement_init_sad = path_to_init_sad + "/displacement_results_dict.pkl"
	
	path_to_strain_results_init_sad = path_to_init_sad + "/local_strain_results_dict.pkl"
	path_to_displacement_init_sad = path_to_init_sad + "/local_displacement_results_dict.pkl"
	
	if os.path.exists(path_to_strain_results_init_sad) and os.path.exists(path_to_displacement_init_sad):
		strains = pickle.load(open(path_to_strain_results_init_sad,'r'))
		disp = pickle.load(open(path_to_displacement_init_sad,'r'))
		local_strains = strains.values()
		local_shear_strain = [strain[1] for strain in local_strains]
		local_vol_strain = [strain[0] for strain in local_strains]
		local_disp = disp.values()
	elif os.path.exists(path_to_all_strain_results_init_sad) and os.path.exists(path_to_all_displacement_init_sad):
		all_strains = pickle.load(open(path_to_all_strain_results_init_sad,'r'))
		all_disp = pickle.load(open(path_to_all_displacement_init_sad,'r'))
		local_strains = [all_strains[k] for k in local_atoms_index if k in all_strains]
		local_shear_strain = [strain[1] for strain in local_strains]
		local_vol_strain = [strain[0] for strain in local_strains]
		local_disp = [all_disp[k] for k in local_atoms_index if k in all_disp]
	else:
		raise Exception("No strain and displacement results found in %s"%path_to_init_sad)
	cluster_ave_shear_strain = np.mean(local_shear_strain)
	cluster_ave_vol_strain = np.mean(local_vol_strain)
	cluster_ave_disp = np.mean(local_disp)
	
	"""
	path_to_sad_fin = path_to_curr_event + "/sad_fin"
	path_to_all_strain_results_sad_fin = path_to_sad_fin + "/strain_results_dict.pkl"
	path_to_all_displacement_sad_fin = path_to_sad_fin + "/displacement_results_dict.pkl"
	
	path_to_strain_results_sad_fin = path_to_sad_fin + "/local_strain_results_dict.pkl"
	path_to_displacement_sad_fin = path_to_sad_fin + "/local_displacement_results_dict.pkl"
	
	if os.path.exists(path_to_strain_results_sad_fin) and os.path.exists(path_to_displacement_sad_fin):
		strains_2 = pickle.load(open(path_to_strain_results_sad_fin,'r'))
		disp_2 = pickle.load(open(path_to_displacement_sad_fin,'r'))
		local_strains_2 = strains_2.values()
		local_shear_strain_2 = [strain[1] for strain in local_strains_2]
		local_vol_strain_2 = [strain[0] for strain in local_strains_2]
		local_disp_2 = disp_2.values()
	elif os.path.exists(path_to_all_strain_results_sad_fin) and os.path.exists(path_to_all_displacement_sad_fin):
		all_strains_2 = pickle.load(open(path_to_all_strain_results_sad_fin,'r'))
		all_disp_2 = pickle.load(open(path_to_all_displacement_sad_fin,'r'))
		local_strains_2 = [all_strains_2[k] for k in local_atoms_index if k in all_strains_2]
		local_shear_strain_2 = [strain[1] for strain in local_strains_2]
		local_vol_strain_2 = [strain[0] for strain in local_strains_2]
		local_disp_2 = [all_disp_2[k] for k in local_atoms_index if k in all_disp_2]
	else:
		raise Exception("No strain and displacement results found in %s"%path_to_sad_fin)
	cluster_ave_shear_strain_2 = np.mean(local_shear_strain_2)
	cluster_ave_vol_strain_2 = np.mean(local_vol_strain_2)
	cluster_ave_disp_2 = np.mean(local_disp_2)
	
	path_to_init_fin = path_to_curr_event + "/init_fin"
	path_to_all_strain_results_init_fin = path_to_init_fin + "/strain_results_dict.pkl"
	path_to_all_displacement_init_fin = path_to_init_fin + "/displacement_results_dict.pkl"
	
	path_to_strain_results_init_fin = path_to_init_fin + "/local_strain_results_dict.pkl"
	path_to_displacement_init_fin = path_to_init_fin + "/local_displacement_results_dict.pkl"
	
	if os.path.exists(path_to_strain_results_init_fin) and os.path.exists(path_to_displacement_init_fin):
		strains_3 = pickle.load(open(path_to_strain_results_init_fin,'r'))
		disp_3 = pickle.load(open(path_to_displacement_init_fin,'r'))
		local_strains_3 = strains_3.values()
		local_shear_strain_3 = [strain[1] for strain in local_strains_3]
		local_vol_strain_3 = [strain[0] for strain in local_strains_3]
		local_disp_3 = disp_3.values()
	elif os.path.exists(path_to_all_strain_results_init_fin) and os.path.exists(path_to_all_displacement_init_fin):
		all_strains_3 = pickle.load(open(path_to_all_strain_results_init_fin,'r'))
		all_disp_3 = pickle.load(open(path_to_all_displacement_init_fin,'r'))
		local_strains_3 = [all_strains_3[k] for k in local_atoms_index if k in all_strains_3]
		local_shear_strain_3 = [strain[1] for strain in local_strains_3]
		local_vol_strain_3 = [strain[0] for strain in local_strains_3]
		local_disp_3 = [all_disp_3[k] for k in local_atoms_index if k in all_disp_3]
	else:
		raise Exception("No strain and displacement results found in %s"%path_to_init_fin)
	cluster_ave_shear_strain_3 = np.mean(local_shear_strain_3)
	cluster_ave_vol_strain_3 = np.mean(local_vol_strain_3)
	cluster_ave_disp_3 = np.mean(local_disp_3)
	"""
	
	return [cluster_ave_shear_strain, cluster_ave_vol_strain, cluster_ave_disp, event]
	
def eng_max_disp(path_to_data_dir, input_param):
	"""
	plot the max displacement vs activiation energy or relaxation energy 
	of all filtered events
	"""
	path_to_act_eng_disp = os.path.join(path_to_data_dir, "act_eng_max_disp.png")
	path_to_relax_eng_disp = os.path.join(path_to_data_dir, "relax_eng_max_disp.png")
	
	path_to_all_act_relax_eng = os.path.join(path_to_data_dir,"act_relax_eng_filtered_events.json")
	if os.path.exists(path_to_all_act_relax_eng):
		saved_results = json.load(open(path_to_all_act_relax_eng, 'r'))
		all_max_disp_A, all_max_disp_B, all_act_eng, all_relax_eng = [],[], [], []
		for result in saved_results:
			event_state = result[0]
			curr_max_disp = event_max_disp(path_to_data_dir, event_state)
			all_max_disp_A.append(curr_max_disp[0])
			all_max_disp_B.append(curr_max_disp[1])
			all_act_eng.append(result[1])
			all_relax_eng.append(result[2])
		plot_2d(path_to_act_eng_disp, all_max_disp_A, all_act_eng, "Max disp /A", "Activation energy /eV")
		plot_2d(path_to_relax_eng_disp, all_max_disp_B, all_relax_eng, "Max disp /A", "Relaxation energy /eV")
	else:
		raise Exception("act_relax_eng_filtered_events.json file does not exist, please do art_data -s SETTINGS --eng --calc first")

def event_max_disp(path_to_data_dir, event_state):
	"""
	get the maximum displacement values for init to sad, sad to fin for
	a single event in data directory
	"""
	path_to_test_dir = os.path.join(path_to_data_dir, event_state[0])
	init, sad, fin = event_state[1][0],event_state[1][1],event_state[1][2]
	path_to_event = path_to_test_dir + "/results/event_" + init + "_" + sad + "_" + fin
	
	max_disp = []
	path_to_init_sad = path_to_event + "/init_sad"
	path_to_displacement_A = path_to_init_sad + "/displacement_results_dict.pkl"
	if os.path.exists(path_to_displacement_A):
		event_disp_A = pickle.load(open(path_to_displacement_A,'r'))
		max_disp_A = max(event_disp_A.values())
		max_disp.append(max_disp_A)
	
	path_to_sad_fin = path_to_event + "/sad_fin"
	path_to_displacement_B = path_to_sad_fin + "/displacement_results_dict.pkl"
	if os.path.exists(path_to_displacement_B):
		event_disp_B = pickle.load(open(path_to_displacement_B,'r'))
		max_disp_B = max(event_disp_B.values())
		max_disp.append(max_disp_B)
	return max_disp
	
def residual_threshold_finder(path_to_data_dir, input_param):
	"""
	this function did a parameter sweep on relative residual threshold
	from 0.01 to 1, the iteration stops when the change of 
	ave_local_atoms[i] - ave_local_atoms[i-1]/relative_threshold[i]-relative_threshold[i-1]
	< -0.2/0.01= -20
	Note: 
	this criteria may be system dependent, need to check for the material system
	
	return:
		ave_local_atoms:
		
	
	"""
	residual_threshold = input_param["residual_threshold"]
	critical_slope = input_param["critical_local_atoms_slope"]
	ave_local_atoms = []
	ave_k = []
	print "all residual_threshold:",residual_threshold
	i=0
	doubt = None
	for x in residual_threshold:
		curr_result = events_local_atoms(path_to_data_dir, input_param, x)
		ave_local_atoms.append(curr_result[0])
		ave_k.append(curr_result[1])
		if i>=1:
			if doubt == True:
				check_criteria = float(ave_local_atoms[i] - ave_local_atoms[i-1])/(residual_threshold[i] - residual_threshold[i-1])
				if check_criteria >= critical_slope:
					print "stopping at relative threshold:", x
					break
				else:
					doubt = False	
			else:
				check_criteria = float(ave_local_atoms[i] - ave_local_atoms[i-1])/(residual_threshold[i] - residual_threshold[i-1])
				if check_criteria >= critical_slope:
					doubt = True
				else:
					doubt = False
		i=i+1
	
	fig, ax1 = plt.subplots()
	ax1.plot(residual_threshold[0:i+1],ave_local_atoms,'rx',markersize=5)
	ax1.set_xlabel('relative residual_threshold')
	ax1.set_ylabel('ave_num_local_atoms')
	
	ax2 = ax1.twinx()
	ax2.plot(residual_threshold[0:i+1], ave_k, 'bo')
	ax2.set_ylabel('ave_slope')
	path_to_image = path_to_data_dir + "/ave_num_k_local_atoms_residual_threshold.png"
	plt.savefig(path_to_image)
	
	print "finally, at relative threshold:%s, average number of locally involved atoms:"%x, ave_local_atoms[-1]
	print "finally, at relative threshold:%s, average slope:"%x, ave_k[-1]
	print "done finding average number of locally involved atoms for current sample at %s"%path_to_data_dir
	return ave_local_atoms[-1], ave_k[-1]


def events_local_atoms_threshold_sweep(path_to_data_dir, input_param):
	# do more on the number of events sweep also and min_sample sweep
	
	# these two parameters are list that containg all parameters need to be sweeped
	#min_samples = input_param["residual_threshold"]
	residual_threshold = input_param["residual_threshold"]
	ave_local_atoms = []
	ave_k = []
	print "residual_threshold:",residual_threshold
	for x in residual_threshold:
		curr_result = events_local_atoms(path_to_data_dir, input_param, x)
		ave_local_atoms.append(curr_result[0])
		ave_k.append(curr_result[1])
	fig, ax1 = plt.subplots()
	ax1.plot(residual_threshold,ave_local_atoms,'rx',markersize=5)
	ax1.set_xlabel('relative residual_threshold')
	ax1.set_ylabel('ave_num_local_atoms')
	
	ax2 = ax1.twinx()
	ax2.plot(residual_threshold, ave_k, 'bo')
	ax2.set_ylabel('ave_slope')
	path_to_image = path_to_data_dir + "/ave_num_k_local_atoms_residual_threshold.png"
	plt.savefig(path_to_image)
	print "done residual threshold parameter sweep for average number of locally involved atoms!"

def all_events_local_atoms_finder(path_to_data_dir, input_param, residual_threshold = 0.5):
	"""
	this function developed correlation model between feature and target using model 
	for all events available in tests with list_of_test_id
	find the outlier atom index and save these local atoms index in a file
	in that event dir
	
	feature: str
		Now only allow "displacement" option
	target: str
		Now allow "shear_strain" option
	Model: str
		Now alow "linear_model" and "LinearSVR" option, which is also adopted when model is None
	"""
	list_of_test_id = input_param["list_of_test_id"]
	model = input_param["model"]
	feature = input_param["feature"]
	target = input_param["target"]
	num_of_proc = input_param["num_of_proc"]
	re_calc = input_param["re_calc"]
	print "current residual_threshold:", residual_threshold
	# perform a function on all events in all tests in list_of_test_id with num_of_proc
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_local_atoms_index(x, path_to_data_dir, model, feature, target, residual_threshold, True, re_calc=re_calc),num_of_proc)
	
	print "done finding all local atoms index for all final selected events in interested tests!"

def all_events_triggered_cluster_atoms_finder(path_to_data_dir, input_param):
	
	list_of_test_id = input_param["list_of_test_id"]
	box_dim = input_param["box_dim"]
	num_of_proc = input_param["num_of_proc"]
	re_calc = input_param["re_calc"]
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_triggered_cluster_atoms_index(x, path_to_data_dir, box_dim, save_results=True, re_calc=re_calc),num_of_proc)
	print "done finding all triggered cluster atoms index for all final selected events in list_of_test_id!"

def single_event_triggered_cluster_atoms_index(event, path_to_data_dir, box_dim, save_results=True, re_calc = False):
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	path_to_curr_result = path_to_test_dir + "/results"
	init,sad,fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	print "path_to_current_event:", path_to_curr_event
	path_to_triggered_atoms_index = path_to_curr_event + "/initial_cluster_atoms_index.json"
	
	if re_calc is False:
		if os.path.exists(path_to_triggered_atoms_index):
			return json.load(open(path_to_triggered_atoms_index,'r'))
	print "re_calculating"
	test_sad_central_atom_id = read_from_art_log_file(path_to_test_dir)
	sad_id = int(sad[3:])
	for curr_sad in test_sad_central_atom_id:
		if curr_sad == sad_id:
			central_atom_id = test_sad_central_atom_id[curr_sad]
			break
	print "current event central atom id:", central_atom_id
	
	path_to_file_ini = path_to_test_dir + '/' + init + ".dump"
	initial_config_data = read_data_from_file(path_to_file_ini)
	central_atom_df = initial_config_data[initial_config_data['item'] == central_atom_id]
	
	cluster_radius = read_cluster_radius_from_bart(path_to_test_dir)
	result_tree = PeriodicCKDTree(box_dim, initial_config_data[['x','y','z']].values)
	result_groups = result_tree.query_ball_point(central_atom_df[['x','y','z']].values, cluster_radius)
	cluster_atoms_df = initial_config_data.iloc[result_groups[0],:]
	cluster_atoms_id = (cluster_atoms_df['item']).tolist()
	
	if save_results is True:
		with open(path_to_triggered_atoms_index, 'w') as f:
			json.dump(cluster_atoms_id,f)
			f.close()
	return cluster_atoms_id
	
	
def events_local_atoms(path_to_data_dir, input_param, residual_threshold = 0.5):
	"""
	this function developed correlation model between feature and target using model 
	for all events available in tests with list_of_test_id
	feature: str
		Now only allow "displacement" option
	target: str
		Now allow "shear_strain" option
	Model: str
		Now alow "linear_model" option, which is also adopted when model is None
	"""
	list_of_test_id = input_param["list_of_test_id"]
	model = input_param["model"]
	feature = input_param["feature"]
	target = input_param["target"]
	num_of_proc = input_param["num_of_proc"]
	print "current residual_threshold:", residual_threshold
	# perform a function on all events in all tests in list_of_test_id with num_of_proc
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_local_atoms(x, path_to_data_dir, model, feature, target, residual_threshold),num_of_proc)
	init_sad_num,sad_fin_num,init_fin_num = [],[],[]
	init_sad_k,sad_fin_k,init_fin_k = [],[],[]
	for event_res in result_list:	
		init_sad_num.append(event_res[0][0])
		sad_fin_num.append(event_res[1][0])
		init_fin_num.append(event_res[2][0])
		init_sad_k.append(event_res[0][1])
		sad_fin_k.append(event_res[1][1])
		init_fin_k.append(event_res[2][1])

	path_to_image_1 = path_to_data_dir + "/num_local_atoms.png"
	plot_histogram_3(path_to_image_1,[init_sad_num,sad_fin_num,init_fin_num])
	ave_num_local_atoms = np.mean(init_sad_num)
	print "the average number of local atoms:", ave_num_local_atoms
	
	path_to_image_2 = path_to_data_dir + "/slope.png"
	plot_histogram_3(path_to_image_2,[init_sad_k,sad_fin_k,init_fin_k])
	#ave_num_local_atoms = sum(init_fin)*1.0/len(init_fin)
	ave_slope = np.mean(init_sad_k)
	print "the average number of slope:", ave_slope
	
	print "done plotting for number of local atoms for all final selected events in interested tests"
	return ave_num_local_atoms, ave_slope
	
def single_event_local_atoms_index(event,path_to_data_dir,model,feature,target,residual_threshold =0.5, save_results=True, re_calc = False):
	"""
	this function developed correlation model between feature and target
	for a single event whose state is saved into event
	
	return:
		local_atom_index: a list
			a list with elements contains a list of local atom index for init_sad
			
	"""
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	path_to_curr_result = path_to_test_dir + "/results"
	init,sad,fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	path_to_init_sad = path_to_curr_event + "/init_sad"
	
	path_to_local_atom_index = path_to_curr_event + "/local_atoms_index.json"
	print "path_to_current_event:", path_to_curr_event
	if re_calc is False:
		if os.path.exists(path_to_local_atom_index):
			return json.load(open(path_to_local_atom_index,'r'))
	print "re_calculating"
	if feature == "displacement" and target == "shear_strain":
		init_sad_X,init_sad_y = get_strain_disp(path_to_init_sad)
	
	init_sad = outlier_detector(init_sad_X,init_sad_y,model,residual_threshold, return_index = True)
	if save_results is True:
		with open(path_to_local_atom_index, 'w') as f:
			json.dump(init_sad,f)
			f.close()
	return init_sad

def single_event_local_atoms(event,path_to_data_dir,model,feature,target,residual_threshold =0.5):
	"""
	this function developed correlation model between feature and target
	for a single event whose state is saved into event
	
	return:
		num_of_local_atom: a list
			num_of_local_atom for init_sad, sad_fin, init_fin
			
	"""
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	path_to_curr_result = path_to_test_dir + "/results"
	init,sad,fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	path_to_init_sad = path_to_curr_event + "/init_sad"
	path_to_sad_fin = path_to_curr_event + "/sad_fin"
	path_to_init_fin = path_to_curr_event + "/init_fin"
	if feature == "displacement" and target == "shear_strain":
		init_sad_X,init_sad_y = get_strain_disp(path_to_init_sad)
		sad_fin_X,sad_fin_y = get_strain_disp(path_to_sad_fin)
		init_fin_X,init_fin_y = get_strain_disp(path_to_init_fin)
	
	init_sad = outlier_detector(init_sad_X,init_sad_y,model,residual_threshold)
	sad_fin = outlier_detector(sad_fin_X,sad_fin_y,model,residual_threshold)
	init_fin = outlier_detector(init_fin_X,init_fin_y,model,residual_threshold)

	return [init_sad, sad_fin, init_fin]

def get_strain_disp(path_to_test_dir):
	path_to_strain_results = path_to_test_dir + "/strain_results_dict.pkl"
	path_to_displacement = path_to_test_dir + "/displacement_results_dict.pkl"
	if os.path.exists(path_to_strain_results) and os.path.exists(path_to_displacement):
		strain = pickle.load(open(path_to_strain_results,'r'))
		displacement = pickle.load(open(path_to_displacement,'r'))
		sample_size = len(strain)
		X =[]
		y =[]
		for key,value in strain.items():
			X.append(displacement[key])
			y.append(value[1])
		X= np.array(X).reshape(sample_size,1)
		y= np.array(y).reshape(sample_size,1)
		return (X,y)
	else:
		raise Exception("strain and displacement data results has not been calculated in %s, can not perform correlation analysis"%path_to_test_dir)
			
	
	
def outlier_detector(feature,target,model=None,residual_threshold = 0.5,return_index = False):
	if model == "linear_model" or model == None:
		return outlier_linear_detector(feature,target,residual_threshold,return_index)
	if model == "LinearSVR":
		# use SVM classification to find inliers, then count outliers
		return outlier_linearSVR_detector(feature,target,residual_threshold,return_index)

def outlier_linear_detector(feature, target, residual_threshold = 0.5, return_index = False):
	"""
	Input argument:
		model: str
			str that specified which sklearn model to use between feature and target
		feature: np.array or pd.dataframe
			the numpy array or pandas.Dataframe that specifies the feature vector or matrix
		target: np.array or pd.dataframe
			the numpy array or pandas.Dataframe that specifies the target value or vector
	Output:
		num_of_outliers:
			the number of outliers that are predicted by the model
		
		r_cut:
		
		
		regression model: sklearn machine leanring class object
			class object of sklearn machine learning model that fits the data
	"""
	
	#base_model = linear_model.LinearRegression()
	#If base_estimator is None, then base_estimator=sklearn.linear_model.LinearRegression() is used 
	# for target values of dtype float
	
	residual_threshold = (np.max(target) -np.min(target))*residual_threshold
	model = linear_model.RANSACRegressor(random_state=1,min_samples=0.1,residual_threshold=residual_threshold)
	model.fit(feature,target)
	inlier_mask = model.inlier_mask_
	outlier_mask = np.logical_not(inlier_mask)
	outlier_index = np.where(outlier_mask)[0]
	num_of_outlier = sum(outlier_mask)
	slope = model.estimator_.coef_[0][0]
	if return_index is False:
		return (num_of_outlier,slope)
	else:
		return outlier_index.tolist()

def outlier_linearSVC_detector(feature,target,residual_threshold, return_index = False):
	"""
	this function detect the outlier by using the SVC with linear kernel
	for this classifer, the y need to be integer array, not float array
	"""
	target = (np.array(target)).flatten()
	#residual_threshold = (np.max(target) -np.min(target))*residual_threshold
	regr = LinearSVC(random_state=1, dual=False, C=residual_threshold,tol=np.mean(target)*0.0001)
	regr.fit(feature, target)
	res_label = regr.predict(feature)
	num_of_outlier = 0
	for label in res_label:
		if label == 0:
			num_of_outlier = num_of_outlier + 1	
	slope = regr.coef_[0]
	return (num_of_outlier, slope)

def outlier_linearSVR_detector(feature, target, residual_threshold, return_index = False):
	"""
	this function detect the outlier by using the LinearSVR with linear kernel
	with the fitted coefficient
	"""
	target = (np.array(target)).flatten()
	residual_threshold = (np.max(target) -np.min(target))*residual_threshold
	regr = LinearSVR(random_state=1, dual=True, epsilon=0.0)
	regr.fit(feature, target)
	
	predict_data = regr.predict(feature)
	i=0
	num_of_outlier = 0
	outlier_index = []
	for x in predict_data:
		delta = x-target[i]
		if abs(delta) > residual_threshold:
			num_of_outlier = num_of_outlier + 1	
			outlier_index.append(i)
		i=i+1
	slope = regr.coef_[0]
	
	if return_index is False:
		return (num_of_outlier, slope)
	else:
		return outlier_index

def fixed_outlier_detector_by_iso_for(feature, outlier_fraction):
	"""
	this function takes a training data X, output the outlier index in X
	with the outlier fraction is outlier_fraction
	"""
	model = IsolationForest(max_samples=100, random_state=1, contamination= outlier_fraction)
	model.fit(feature)
	y_predict = model.predict(feature)
	outliers = []
	i=0
	for y in y_predict:
		if y == -1:
			outliers.append(i)
		i=i+1
	return outliers

def fixed_outlier_detector_by_LOF(feature, outlier_fraction):
	"""
	this function takes a training data X, output the outlier index in X
	with the outlier fraction is outlier_fraction
	"""
	model = NearestNeighbors(contamination= outlier_fraction)
	model.fit(feature)
	y_predict = model.predict(feature)
	outliers = []
	i=0
	for y in y_predict:
		if y == -1:
			outliers.append(i)
		i=i+1
	return outliers

def feature_to_X(feature):
	"""
	this function converts the original data into X with shape (n_sample,n_feature)
	that will be used by sklearn class object
	"""

def feature_to_y(target):
	"""
	this function converts the original target data into y with shape (n_sample,n_target)
	that will be used by sklearn class object
	"""
	


def basin_config_distance_activation_energy(path_to_test_dir, event_state):
	"""
	each test is corresponding to the triggering of a single atom, which search the
	potential energy trajectory around that basin of that atom located
	
	event_state is a list of ordered states (corresponding to their appearance order), 
	each event state contains the strings of all states, init_state, sad_state, fin_state, 
	with intermediate states that can be obtained from ART output
	
	Input:
		
		event_state: a list
			a list of strs with each str being an state, e,g min1000
	
	Returns:
	
		config_distance: a list
			a list of distances from current state to the initial state
		
		eng_barrier: a list
			a list of energy barriers from current state to the initial state
	"""
	i=0
	config_distance, eng_barrier = []
	for state in event_state:
		
		if i == 0:
			init = state
			path_to_init_file = path_to_test_dir + '/' + init + ".dump"
			init_config = Configuration(path_to_init_file)
			config_distance.append(0)
			eng_barrier.append(0)
		else:
			path_to_state_file = path_to_test_dir + '/' + state + ".dump"
			
			state_config = Configuration(path_to_state_file)
			
			state_act_eng = state_energy_barrier(path_to_test_dir, init, state)
			
			state_distance = Configuration.distance_pbc(state_config, init_config)
			
			config_distance.append(state_distance)
			
			eng_barrier.append(state_act_eng)
	return (config_distance, eng_barrier)

