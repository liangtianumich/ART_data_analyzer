import os, json, pickle
import numpy as np
import pandas as pd
import time
import random
import multiprocessing as mp
import argparse
from util import prompt_yes_no, operation_on_events, data_dir_to_test_dir,get_list_of_atoms_from_atom_list, event_act_relax_energy, event_ave_strain_displacement, event_voronoi_volume, event_voronoi_class
from data_reader import read_data_from_file


def generate_correlation_table_mp(path_to_data_dir, input_param):
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	atom_list = input_param["atom_list"]
	print "if atom_list is local or pn or pn_non_affine, voronoi volume calculation only act on atom_list from init to sad"
	print "confirm if voronoi_index_results.json is corresponding to the atom_list you just specified:", atom_list
	if not prompt_yes_no():
		raise Exception("quitting, re_calc the voronoi indexes for your specified atom_list by --voro --calc --re_calc --local if atom_list is local")
	
	operation = lambda x: single_event_data_extractor(x, path_to_data_dir,atom_list)
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	convert_to_csv(path_to_data_dir, result_list)
	print "All done!"
	

def single_event_data_extractor(event_state, path_to_data_dir, atom_list):
	# go to event directory
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	path_to_curr_result = path_to_test_dir + "/results"	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
		
	print "extracting atom_list:"
	path_to_file_ini = path_to_test_dir + '/' + init + ".dump"
	initial_config_data = read_data_from_file(path_to_file_ini)
	(atom_list_init_sad, atom_list_sad_fin, atom_list_init_fin) = get_list_of_atoms_from_atom_list(path_to_curr_event, initial_config_data, atom_list)	
	atom_num_init_sad, atom_num_sad_fin, atom_num_init_fin = len(atom_list_init_sad),len(atom_list_sad_fin),len(atom_list_init_fin)
	
	print "extracting energy data:"
	event, act_eng, relax_eng = event_act_relax_energy(event_state, path_to_data_dir)
	init_sad_eng, sad_fin_eng, init_fin_eng = act_eng, -relax_eng, act_eng - relax_eng
	
	print "extracting strain and displacement data:"
	# strain = {"init_sad":[], "sad_fin":[], "init_fin:[]"}
	# displacement = {"init_sad":[], "sad_fin":[], "init_fin:[]"}
	# atom_list= {"init_sad":[], "sad_fin":[], "init_fin:[]"}
	atom_list= {"init_sad":atom_list_init_sad, "sad_fin":atom_list_sad_fin, "init_fin":atom_list_init_fin}
	ave_vol_strain, ave_shear_strain, ave_disp = event_ave_strain_displacement(event_state, path_to_data_dir, atom_list)

	print "extracting voronoi volume data:"
	init_vol, sad_vol, fin_vol = event_voronoi_volume(event_state, path_to_data_dir)
	init_sad_vol, sad_fin_vol, init_fin_vol = sad_vol - init_vol, fin_vol - sad_vol, fin_vol - init_vol
	
	event_voro_res = event_voronoi_class(event_state, path_to_data_dir)
	
	event_init_sad = pd.DataFrame(dtype=object)
	event_init_sad.at[0,"event_id"] = str(event_state)
	event_init_sad.at[0,"atom_list"] = str(atom_list_init_sad)
	event_init_sad.at[0,"num_of_atoms"] = str(atom_num_init_sad)
	event_init_sad.at[0,"ave_vol_strain"] = str(ave_vol_strain["init_sad"])
	event_init_sad.at[0,"ave_shear_strain"] = str(ave_shear_strain["init_sad"])
	event_init_sad.at[0,"ave_disp"] = str(ave_disp["init_sad"])
	event_init_sad.at[0,"start_ICO_frac"] = str(event_voro_res["init"][0])
	event_init_sad.at[0,"start_ICO_like_frac"] = str(event_voro_res["init"][1])
	event_init_sad.at[0,"start_GUM_frac"] = str(event_voro_res["init"][2])
	event_init_sad.at[0,"end_ICO_frac"] = str(event_voro_res["sad"][0])
	event_init_sad.at[0,"end_ICO_like_frac"] = str(event_voro_res["sad"][1])
	event_init_sad.at[0,"end_GUM_frac"] = str(event_voro_res["sad"][2])
	event_init_sad.at[0,"vol_diff"] = str(init_sad_vol)
	event_init_sad.at[0,"eng_diff"] = str(init_sad_eng)
	
	event_sad_fin = pd.DataFrame(dtype=object)
	event_sad_fin.at[0,"event_id"] = str(event_state)
	event_sad_fin.at[0,"atom_list"] = str(atom_list_sad_fin)
	event_sad_fin.at[0,"num_of_atoms"] = str(atom_num_sad_fin)
	event_sad_fin.at[0,"ave_vol_strain"] = str(ave_vol_strain["sad_fin"])
	event_sad_fin.at[0,"ave_shear_strain"] = str(ave_shear_strain["sad_fin"])
	event_sad_fin.at[0,"ave_disp"] = str(ave_disp["sad_fin"])
	event_sad_fin.at[0,"start_ICO_frac"] = str(event_voro_res["sad"][0])
	event_sad_fin.at[0,"start_ICO_like_frac"] = str(event_voro_res["sad"][1])
	event_sad_fin.at[0,"start_GUM_frac"] = str(event_voro_res["sad"][2])
	event_sad_fin.at[0,"end_ICO_frac"] = str(event_voro_res["fin"][0])
	event_sad_fin.at[0,"end_ICO_like_frac"] = str(event_voro_res["fin"][1])
	event_sad_fin.at[0,"end_GUM_frac"] = str(event_voro_res["fin"][2])
	event_sad_fin.at[0,"vol_diff"] = str(sad_fin_vol)
	event_sad_fin.at[0,"eng_diff"] = str(sad_fin_eng)
	
	event_init_fin = pd.DataFrame(dtype=object)
	event_init_fin.at[0,"event_id"] = str(event_state)
	event_init_fin.at[0,"atom_list"] = str(atom_list_init_fin)
	event_init_fin.at[0,"num_of_atoms"] = str(atom_num_init_fin)
	event_init_fin.at[0,"ave_vol_strain"] = str(ave_vol_strain["init_fin"])
	event_init_fin.at[0,"ave_shear_strain"] = str(ave_shear_strain["init_fin"])
	event_init_fin.at[0,"ave_disp"] = str(ave_disp["init_fin"])
	event_init_fin.at[0,"start_ICO_frac"] = str(event_voro_res["init"][0])
	event_init_fin.at[0,"start_ICO_like_frac"] = str(event_voro_res["init"][1])
	event_init_fin.at[0,"start_GUM_frac"] = str(event_voro_res["init"][2])
	event_init_fin.at[0,"end_ICO_frac"] = str(event_voro_res["fin"][0])
	event_init_fin.at[0,"end_ICO_like_frac"] = str(event_voro_res["fin"][1])
	event_init_fin.at[0,"end_GUM_frac"] = str(event_voro_res["fin"][2])
	event_init_fin.at[0,"vol_diff"] = str(init_fin_vol)
	event_init_fin.at[0,"eng_diff"] = str(init_fin_eng)
		
	#event_init_sad = = pd.DataFrame({"event_id":event_state,"atom_list":atom_list_init_sad,"num_of_atoms":atom_num_init_sad,"ave_vol_strain":ave_vol_strain["init_sad"],"ave_shear_strain":ave_shear_strain["init_sad"],"ave_disp":ave_disp["init_sad"],"start_ICO_frac":event_voro_res["init"][0],"start_ICO_like_frac":event_voro_res["init"][1],"start_GUM_frac":event_voro_res["init"][2],"end_ICO_frac":event_voro_res["sad"][0],"end_ICO_like_frac":event_voro_res["sad"][1],"end_GUM_frac":event_voro_res["sad"][2], "vol_diff":init_sad_vol,"eng_diff":init_sad_eng})
	#event_sad_fin = pd.DataFrame({"event_id":event_state,"atom_list":atom_list_sad_fin,"num_of_atoms":atom_num_sad_fin,"ave_vol_strain":ave_vol_strain["sad_fin"],"ave_shear_strain":ave_shear_strain["sad_fin"],"ave_disp":ave_disp["sad_fin"],"start_ICO_frac":event_voro_res["sad"][0],"start_ICO_like_frac":event_voro_res["sad"][1],"start_GUM_frac":event_voro_res["sad"][2],"end_ICO_frac":event_voro_res["fin"][0],"end_ICO_like_frac":event_voro_res["fin"][1],"end_GUM_frac":event_voro_res["fin"][2],"vol_diff": sad_fin_vol,"eng_diff":sad_fin_eng})
	#event_init_fin = pd.DataFrame({"event_id":event_state,"atom_list":atom_list_init_fin,"num_of_atoms":atom_num_init_fin,"ave_vol_strain":ave_vol_strain["init_fin"],"ave_shear_strain":ave_shear_strain["init_fin"],"ave_disp":ave_disp["init_fin"],"start_ICO_frac":event_voro_res["init"][0],"start_ICO_like_frac":event_voro_res["init"][1],"start_GUM_frac":event_voro_res["init"][2],"end_ICO_frac":event_voro_res["fin"][0],"end_ICO_like_frac":event_voro_res["fin"][1],"end_GUM_frac":event_voro_res["fin"][2],"vol_diff": init_fin_vol,"eng_diff":init_fin_eng})
	
	return (event_init_sad,event_sad_fin,event_init_fin)
	
def convert_to_csv(path_to_data_dir, result_list):
	init_sad_pd, sad_fin_pd, init_fin_pd = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
	for result in result_list:
		init_sad_pd = init_sad_pd.append(result[0], ignore_index=True)
		sad_fin_pd = sad_fin_pd.append(result[1], ignore_index=True)
		init_fin_pd = init_fin_pd.append(result[2], ignore_index=True)
	
	path_to_init_sad = os.path.join(path_to_data_dir,"init_sad_data.csv")
	path_to_sad_fin = os.path.join(path_to_data_dir,"sad_fin_data.csv")
	path_to_init_fin = os.path.join(path_to_data_dir,"init_fin_data.csv")
	
	init_sad_pd.to_csv(path_or_buf=path_to_init_sad)
	sad_fin_pd.to_csv(path_or_buf=path_to_sad_fin)
	init_fin_pd.to_csv(path_or_buf=path_to_init_fin)
		
		
		
