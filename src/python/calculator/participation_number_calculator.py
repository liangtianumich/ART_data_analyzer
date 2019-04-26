import os
import math
import numpy as np
import pickle
import json
from util import operation_on_events, data_dir_to_test_dir

def pn_calculator_run_all_tests_mp(path_to_data_dir, input_param, save_results = True, is_non_affine = False):
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	re_calc = input_param["re_calc"]
	
	operation = lambda x: single_event_pn_calculator(x, path_to_data_dir, re_calc = re_calc, is_non_affine = is_non_affine)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	print ("done pn calculations for all interested tests!")
	return result_list

def single_event_pn_calculator(event_state, path_to_data_dir, save_results = True, re_calc = False, is_non_affine = False):
	# go to event directory
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	path_to_curr_result = path_to_test_dir + "/results"	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	if not os.path.exists(path_to_curr_event):
		os.makedirs(path_to_curr_event)
	# check if results has already been saved
	path_pn_number_results = os.path.join(path_to_curr_event, "pn_number.json")
	path_pn_index_results = os.path.join(path_to_curr_event, "pn_index.json")
		
	if re_calc is False:
		if os.path.exists(path_pn_number_results) and os.path.exists(path_pn_index_results):
			return (json.load(open(path_pn_number_results,"r")), json.load(open(path_pn_index_results,"r")))
		else:
			print "begin calculating pn number and index"
	else:
		print "re_calculating pn number and index"
	
	# begin working on each process of this event
	path_to_init_sad = path_to_curr_event + "/init_sad"
	path_to_sad_fin = path_to_curr_event + "/sad_fin"
	path_to_init_fin = path_to_curr_event + "/init_fin"
	
	if is_non_affine is False:
		# read displacement file for each event
		path_to_all_displacement_init_sad = path_to_init_sad + "/displacement_results_dict.pkl"
		path_to_all_displacement_sad_fin = path_to_sad_fin + "/displacement_results_dict.pkl"
		path_to_all_displacement_init_fin = path_to_init_fin + "/displacement_results_dict.pkl"
		if os.path.exists(path_to_all_displacement_init_sad) and os.path.exists(path_to_all_displacement_sad_fin) and os.path.exists(path_to_all_displacement_init_fin):
			all_disp_init_sad = pickle.load(open(path_to_all_displacement_init_sad,'r'))
			all_disp_sad_fin = pickle.load(open(path_to_all_displacement_sad_fin,'r'))
			all_disp_init_fin = pickle.load(open(path_to_all_displacement_init_fin,'r'))
		else:
			raise Exception("displacement results have not been found in current event:%s"%path_to_curr_event)
		# perform the calculations according to the formulas proposed in 
		# "Local structural excitations in model glasses" Swayamjyoti et al PRB, 2014
		pn_init_sad,pn_index_init_sad = get_pn_number_index(all_disp_init_sad)
		pn_sad_fin,pn_index_sad_fin = get_pn_number_index(all_disp_sad_fin)
		pn_init_fin,pn_index_init_fin = get_pn_number_index(all_disp_init_fin)
	elif is_non_affine is True:
		path_to_all_shear_init_sad = path_to_init_sad + "/strain_results_dict.pkl"
		path_to_all_shear_sad_fin = path_to_sad_fin + "/strain_results_dict.pkl"
		path_to_all_shear_init_fin = path_to_init_fin + "/strain_results_dict.pkl"
		if os.path.exists(path_to_all_shear_init_sad) and os.path.exists(path_to_all_shear_sad_fin) and os.path.exists(path_to_all_shear_init_fin):
			all_strains_init_sad = pickle.load(open(path_to_all_shear_init_sad,'r'))
			all_strains_sad_fin = pickle.load(open(path_to_all_shear_sad_fin,'r'))
			all_strains_init_fin = pickle.load(open(path_to_all_shear_init_fin,'r'))
			
			all_shear_init_sad = extract_shear_dict(all_strains_init_sad)
			all_shear_sad_fin = extract_shear_dict(all_strains_sad_fin)
			all_shear_init_fin = extract_shear_dict(all_strains_init_fin)
			
		else:
			raise Exception("strain results have not been found in current event:%s"%path_to_curr_event)
		# perform the calculations according to the formulas proposed in 
		# "Local structural excitations in model glasses" Swayamjyoti et al PRB, 2014
		pn_init_sad,pn_index_init_sad = get_pn_number_index(all_shear_init_sad,True)
		pn_sad_fin,pn_index_sad_fin = get_pn_number_index(all_shear_sad_fin,True)
		pn_init_fin,pn_index_init_fin = get_pn_number_index(all_shear_init_fin,True)
	
	
	pn_number = {"init_sad": pn_init_sad, "sad_fin":pn_sad_fin, "init_fin": pn_init_fin}
	pn_index = {"init_sad": pn_index_init_sad, "sad_fin":pn_index_sad_fin, "init_fin": pn_index_init_fin}
	# save the calculated results into files
	if save_results is True:
		print ("begin saving pn results into json files")
		with open(path_pn_number_results, 'w') as f:
			json.dump(pn_number,f)
			f.close()
		with open(path_pn_index_results, 'w') as f:
			json.dump(pn_index,f)
			f.close()
		print "pn results saved into two json files!"
	return (pn_number,pn_index)
	
def get_pn_number_index(all_disp_dict,is_non_affine=False):
	if is_non_affine is True:
		for key, value in all_disp_dict.items():
			all_disp_dict[key] = abs(value)

	disp_list = all_disp_dict.values()
	weights = []
	for disp in disp_list:
		weights.append(disp ** 4)
	total_weights = sum(weights)
	norm_weights = [x*1.0/total_weights for x in weights]
	
	norm_sq = []
	for y in norm_weights:
		norm_sq.append(y ** 2)
	
	pn_number = math.ceil(1.0/sum(norm_sq))
	
	max_indexes = np.argsort(all_disp_dict.values())[-int(pn_number):]
	
	pn_indexes = [all_disp_dict.keys()[i] for i in max_indexes]
	return (pn_number,pn_indexes)
			
def extract_shear_dict(all_strains):
	for key, value in all_strains.items():
		all_strains[key] = value[1]
	return all_strains

	

	

		
