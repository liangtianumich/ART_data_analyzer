"""
this module implements correlation model
"""
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from sklearn import linear_model
from visualizer.strain_visualizer import plot_histogram_3
from util import operation_on_events, Configuration, state_energy_barrier

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
	#plt.figure()
	#path_to_image = path_to_data_dir + "/ave_num_local_atoms_residual_threshold.png"
	#plt.plot(residual_threshold,ave_local_atoms,'rx',markersize=5)
	#plt.xlabel('residual_threshold',fontsize=20)
	#plt.ylabel('ave_num_local_atoms',fontsize=20)
	#plt.savefig(path_to_image)
	#plt.close()
	print "done residual threshold parameter sweep for average number of locally involved atoms!"

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
	
	path_to_image_2 = path_to_data_dir + "/slope.png"
	plot_histogram_3(path_to_image_2,[init_sad_k,sad_fin_k,init_fin_k])
	#ave_num_local_atoms = sum(init_fin)*1.0/len(init_fin)
	ave_num_local_atoms = np.mean(init_sad_num)
	ave_slope = np.mean(init_sad_k)
	print "the average number of local atoms:", ave_num_local_atoms
	print "the average number of slope:", ave_slope
	print "done plotting for number of local atoms for all final selected events in interested tests"
	return ave_num_local_atoms, ave_slope

def single_event_local_atoms(event,path_to_data_dir,model,feature,target,residual_threshold =0.5):
	"""
	this function developed correlation model between feature and target
	for a single event whose state is saved into event
	
	return:
		num_of_local_atom: a list
			num_of_local_atom for init_sad, sad_fin, init_fin
			
	"""
	path_to_curr_result = path_to_data_dir + event[0] + "/results"
	init,sad,fin = event[1][0], event[1][1], event[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	path_to_init_sad = path_to_curr_event + "/init_sad"
	path_to_sad_fin = path_to_curr_event + "/sad_fin"
	path_to_init_fin = path_to_curr_event + "/init_fin"
	if feature == "displacement" and target == "shear_strain":
		init_sad_X,init_sad_y = get_strain_disp(path_to_init_sad)
		sad_fin_X,sad_fin_y = get_strain_disp(path_to_sad_fin)
		init_fin_X,init_fin_y = get_strain_disp(path_to_init_fin)
	
	#path_to_init_sad_strain_results = path_to_init_sad + "/strain_results_dict.pkl"
	#path_to_init_sad_displacement = path_to_init_sad + "/displacement_results_dict.pkl"
	
	#path_to_sad_fin_strain_results = path_to_sad_fin + "/strain_results_dict.pkl"
	#path_to_sad_fin_displacement = path_to_sad_fin + "/displacement_results_dict.pkl"
	
	#path_to_init_fin_strain_results = path_to_init_fin + "/strain_results_dict.pkl"
	#path_to_init_fin_displacement = path_to_init_fin + "/displacement_results_dict.pkl"
	init_sad = outlier_detector(init_sad_X,init_sad_y,model,residual_threshold)
	sad_fin = outlier_detector(sad_fin_X,sad_fin_y,model,residual_threshold)
	init_fin = outlier_detector(init_fin_X,init_fin_y,model,residual_threshold)

	return [init_sad,sad_fin,init_fin]
	
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
			
	
	
def outlier_detector(feature,target,model=None,residual_threshold = 0.5):
	if model == "linear_model" or model == None:
		return outlier_linear_detector(feature,target,residual_threshold)
	#if model == "SVM":
		# use SVM classification to find inliers, then count outliers
		#

def outlier_linear_detector(feature, target, residual_threshold = 0.5):
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
	num_of_outlier = sum(outlier_mask)
	slope = model.estimator_.coef_[0][0]
	return (num_of_outlier,slope)

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

