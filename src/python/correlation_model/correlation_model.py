"""
this module implements correlation model
"""
import os
import pickle
import numpy as np
import multiprocessing as mp
from sklearn import linear_model
from visualizer.strain_visualizer import plot_histogram_3
from util import operation_on_events

def events_correlation_model(path_to_data_dir, input_param):
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
	
	# perform a function on all events in all tests in list_of_test_id with num_of_proc
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, lambda x: single_event_correlation(x,path_to_data_dir,model,feature,target),num_of_proc)
	init_sad,sad_fin,init_fin = [],[],[]
	for event_res in result_list:
		init_sad.append(event_res[0])
		sad_fin.append(event_res[1])
		init_fin.append(event_res[2])
	path_to_image = path_to_data_dir + "/num_local_atoms.png"
	plot_histogram_3(path_to_image,[init_sad,sad_fin,init_fin])
	print "the average number of local atoms:", sum(init_fin)*1.0/len(init_fin)
	print "done plotting for number of local atoms for all final selected events in interested tests"

def single_event_correlation(event,path_to_data_dir,model,feature,target):
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
	init_sad_num = outlier_correlation_model(init_sad_X,init_sad_y,model)
	sad_fin_num = outlier_correlation_model(sad_fin_X,sad_fin_y,model)
	init_fin_num = outlier_correlation_model(init_fin_X,init_fin_y,model)

	return [init_sad_num,sad_fin_num,init_fin_num]
	
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
			
	
	
def outlier_correlation_model(feature,target,model=None):
	if model == "linear_model" or model == None:
		return outlier_linear_model(feature,target)
	#if model == "SVM":
		# use SVM classification to find inliers, then count outliers
		#

def outlier_linear_model(feature,target,model=None):
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
	if model is None:
		model = linear_model.RANSACRegressor(random_state=1)
	model.fit(feature,target)
	inlier_mask = model.inlier_mask_
	outlier_mask = np.logical_not(inlier_mask)
	num_of_outlier = sum(outlier_mask)
	return num_of_outlier

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
	

def outlier_detector(model,X,y):
	"""
	Input arguments:
		model: scikit-learn class
			choice of scikit-learn class object
		X : numpy array or sparse matrix of shape [n_samples,n_features] Training data
		y: numpy array of shape [n_samples,n_targets] Target values Will be cast to dtype if necessary
	return:
		number of outlier points:
		
		r_cut:
		
		index of data outlier points:
	"""
	


def PEL_config_distance_activiation_energy(path_to_test_dir, event_state):
	"""
	event_state is a list containing the strings of all states,
	init_state, sad_state, fin_state, with intermediate states
	
	"""
	# trajectory
	
	#Configuration(path_to_init_file)
	
	#Configuration(path_to_sad_file)
	
	#Cinfiguration(path_to_fin_file)

