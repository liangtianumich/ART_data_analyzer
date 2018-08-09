"""
this event energy calculator module extract various energy related quantities from
each event
"""
from util import event_act_relax_energy, get_list_of_final_filtered_events_str
import multiprocessing as mp
import pandas as pd
import numpy as np
import os, json
from functools import partial
from scipy.stats import ttest_ind, ttest_rel
from visualizer.event_energy_visualizer import plot_act_relax_histogram
from sklearn.model_selection import KFold

def energy_calculator_run_all_tests_mp(path_to_data_dir, input_param, save_results = True, re_calc = False):
	"""
	this function extract the activation energy and relaxation energy
	for all final filtered events in all tests in list_of_test_id, summary of their statistics
	, for the purpose of checking the convergence of their distribution
	"""
	list_of_test_id = input_param['list_of_test_id']
	
	if 'num_of_proc' in input_param and 're_calc' in input_param:
		num_of_proc = input_param['num_of_proc']
		re_calc = input_param["re_calc"]
	
	path_to_all_act_relax_eng = path_to_data_dir + "act_relax_eng_filtered_events.json"
	
	if re_calc is False:
		if os.path.exists(path_to_all_act_relax_eng):
			saved_results = json.load(open(path_to_all_act_relax_eng, 'r'))
			all_act_eng, all_relax_eng = [], []
			for result in saved_results:
				all_act_eng.append(result[1])
				all_relax_eng.append(result[2])
			return [all_act_eng, all_relax_eng]
	
	list_event_str = get_list_of_final_filtered_events_str(path_to_data_dir)
	pool = mp.Pool(processes = num_of_proc)
	result = pool.map(partial(event_act_relax_energy,path_to_data_dir = path_to_data_dir), list_event_str)
	
	all_act_eng, all_relax_eng = [], []
	saved_data = []
	for event_result in result:
		(event, act_eng, relax_eng) = event_result
		saved_data.append([event,act_eng,relax_eng])
		all_act_eng.append(act_eng)
		all_relax_eng.append(relax_eng)
	
	# get statistics
	energy_results = pd.DataFrame({"act_eng":np.array(all_act_eng),"relax_eng": np.array(all_relax_eng)})
	energy_results.describe()
	
	# plot the energy distribution
	path_to_eng_plot = path_to_data_dir + "act_relax_eng_filtered_events.png"
	
	plot_act_relax_histogram(path_to_eng_plot, [all_act_eng, all_relax_eng])
	print "done finding all activation energy and relaxation energy for all filtered events in list_of_test_id"
	
	if save_results is True:
		json.dump(saved_data, open(path_to_all_act_relax_eng,'w'))

	return [all_act_eng, all_relax_eng]


def eng_convergence_ttest_ind(path_to_data_dir_1, path_to_data_dir_2, equal_var = True):
	
	"""
	this function takes the act_eng and relax_eng data from two independent samples in 
	two ART data directories with different input parameters, perform the two independent sample t tests
	to check if the mean of samples are identical to confirm if the act and relax energy distributions converge
	
	return True if two data has identical average
	return False if two data does not have identical average
	"""
	path_to_eng_1 = path_to_data_dir_1 + "act_relax_eng_filtered_events.json"
	path_to_eng_2 = path_to_data_dir_2 + "act_relax_eng_filtered_events.json"
	if os.path.exists(path_to_eng_1) and os.path.exists(path_to_eng_2):
		eng_data_1 = json.load(open(path_to_eng_1, 'r'))
		eng_data_2 = json.load(open(path_to_eng_2, 'r'))
	else:
		raise Exception("at least one of the data directories do not have the energy data file")
	act_eng_1, relax_eng_1 = [], []
	for event in eng_data_1:
		act_eng_1.append(event[1])
		relax_eng_1.append(event[2])
	
	act_eng_2, relax_eng_2 = [], []
	for event in eng_data_2:
		act_eng_2.append(event[1])
		relax_eng_2.append(event[2])
	
	t_1, prob_1 = ttest_ind(np.array(act_eng_1), np.array(act_eng_2), equal_var = equal_var)
	
	t_2, prob_2 = ttest_ind(np.array(relax_eng_1), np.array(relax_eng_2), equal_var = equal_var)
	
	if prob_1 < 0.05 and prob_2 < 0.05:
		return True
	else:
		return False

def eng_convergence_ttest_rel(path_to_data_dir_1, path_to_data_dir_2):
	
	"""
	this function takes the act_eng and relax_eng data from two related samples 
	(two different sampled data subsets with the same length) from the same LAMMPS sample with 
	the same ART input parameters, perform the related t tests 
	to check if the mean of samples are identical to confirm if the act 
	and relax energy distributions converge
	
	return True if two subset data has identical average
	return False if two subset data does not have identical average
	"""
	path_to_eng_1 = path_to_data_dir_1 + "act_relax_eng_filtered_events.json"
	path_to_eng_2 = path_to_data_dir_2 + "act_relax_eng_filtered_events.json"
	if os.path.exists(path_to_eng_1) and os.path.exists(path_to_eng_2):
		eng_data_1 = json.load(open(path_to_eng_1, 'r'))
		eng_data_2 = json.load(open(path_to_eng_2, 'r'))
	else:
		raise Exception("at least one of the data directories do not have the energy data file")
	act_eng_1, relax_eng_1 = [], []
	for event in eng_data_1:
		act_eng_1.append(event[1])
		relax_eng_1.append(event[2])
	
	act_eng_2, relax_eng_2 = [], []
	for event in eng_data_2:
		act_eng_2.append(event[1])
		relax_eng_2.append(event[2])
	
	t_1, prob_1 = ttest_rel(np.array(act_eng_1), np.array(act_eng_2))
	
	t_2, prob_2 = ttest_rel(np.array(relax_eng_1), np.array(relax_eng_2))
	
	if prob_1 < 0.05 and prob_2 < 0.05:
		return True
	else:
		return False

def eng_k_fold_ttest(path_to_data_dir, k=2, option='ind', n=1):
	"""
	A more rigorous way to check convergence
	this function randomly divides the data in path_to_data_dir into k_folds in n times,
	for each time, do t test on each fold-pair of the k folder,
	if all t tests satisfied the criteria, then we can confidently say that
	the amount of eng data is convergent
	"""
	print "current t test mode is %s"%option
	path_to_eng = path_to_data_dir + "act_relax_eng_filtered_events.json"
	
	if os.path.exists(path_to_eng):
		eng_data = json.load(open(path_to_eng, 'r'))
	else:
		raise Exception("the data directories %s do not have the energy data file"%path_to_data_dir)
	
	act_eng, relax_eng = [], []
	for event in eng_data:
		act_eng.append(event[1])
		relax_eng.append(event[2])
	
	# or we can kfold act_eng and relax_eng independently by KFold_n_times(X,k,n)
	kfolds_n_act,kfolds_n_relax = KFold_xy_n_times(act_eng,relax_eng,k,n)
	
	for i in range(n):
		kfolds_act, kfolds_relax = kfolds_n_act[i],kfolds_n_relax[i]
		print "k fold %s times"%(i+1)
		is_act, is_relax = ttest_kfold(kfolds_act, option=option),ttest_kfold(kfolds_relax, option=option)
		if is_act == False or is_relax == False:
			return False
	return True
			

def ttest_kfold(kfolds, option='ind'):
	if option == 'ind':
		return ttest_ind_kfold(kfolds)
	elif option == 'rel':
		return ttest_rel_kfold(kfolds)
			
def ttest_ind_kfold(kfolds):
	"""
	this function perform t test on all possible fold pairs, return True
	if all fold pairs are convergent, return False otherwise
	kfolds: list or numpy array
	"""
	k = len(kfolds)
	for i in range(k):
		for j in range(i+1,k):
			t_1, prob_1 = ttest_ind(np.array(kfolds[i]), np.array(kfolds[j]))
			if prob_1 > 0.05:
				return False
	return True

def ttest_rel_kfold(kfolds):
	"""
	this function perform t test on all possible fold pairs, return True
	if all fold pairs are convergent, return False otherwise
	kfolds: list or numpy array
	"""
	k = len(kfolds)
	for i in range(k):
		for j in range(i+1,k):
			t_1, prob_1 = ttest_rel(np.array(kfolds[i]), np.array(kfolds[j]))
			if prob_1 > 0.05:
				return False
	return True
	
def KFold_n_times(X,k,n):
	"""
	this function returns a list, in which each element is a list that 
	contains the numpy array of each fold when splitting X data into k folds
	"""
	X = np.array(X)
	kfolds_n = []
	for i in range(n):
		kf = KFold(n_splits=k, shuffle=True)
		kfolds = []
		for train_index, test_index in kf.split(X):
			kfolds.append(X[test_index])
		kfolds_n.append(kfolds)
	return kfolds_n

def KFold_xy_n_times(X,y,k,n):
	"""
	this function returns a tuple, for each element of the tuple, is a list that 
	contains the numpy array of each fold when splitting X or y data into k folds
	"""
	X = np.array(X)
	y = np.array(y)
	kfolds_n_X, kfolds_n_y = [], []
	for i in range(n):
		kf = KFold(n_splits=k, shuffle=True)
		kfolds_x, kfolds_y = [],[]
		for train_index, test_index in kf.split(X):
			kfolds_x.append(X[test_index])
			kfolds_y.append(y[test_index])
		kfolds_n_X.append(kfolds_x)
		kfolds_n_y.append(kfolds_y)
	return (kfolds_n_X,kfolds_n_y)
