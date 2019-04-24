"""
this module runs ART simulations by set up the ART input files for each MD dump file sample, 
create a run_dir associated with the central_atom id, 
copy all modified input files into each run_dir, run ./bart.sh under each run_dir
by subprocess module in parallel
"""
import subprocess
import os
import re
import json
import shutil
import copy
import multiprocessing as mp
from data_reader import *
from util import data_dir_to_test_dir, prompt_yes_no, update_central_atom_list


def run_art_mp(path_to_data_dir, input_param, pbs=False):
	"""
	this function will run art by bart.sh inside the path_to_data_dir for each path_to_run_dir
	in a parallel fashion
	"""
	
	num_of_proc = input_param["num_of_proc"]
	central_atom_list = input_param['central_atom_list']
	
	pool = mp.Pool(processes=num_of_proc)
		
	list_of_run_dir = []
	for central_atom in central_atom_list:
		path_to_run_dir = os.path.join(path_to_data_dir, str(central_atom))
		list_of_run_dir.append(path_to_run_dir)
	pool.map(run_bart,list_of_run_dir)
	print "done running ART!"

def run_bart(path_to_run_dir):
	# run bart.sh inside each dir
	os.chdir(path_to_run_dir)
	subprocess.check_call("chmod +x mod_bart.sh; ./mod_bart.sh", shell=True)
	print "done runing art for %s"%path_to_run_dir
	

def set_up_input_files(path_to_data_dir, input_param):
	"""
	this function will set up the input files for each MD dump file sample,
	, create a run_dir associated with the central_atom id,
	copy all modified input files into each run_dir,
	run bart.sh under each run_dir
	
	directory structure: 
	
	# path_to_input_files are saved in the $ART_INPUT
	# path_to_data_dir = path_to_root_dir + sample_name
	# path_to_run_dir = path_to_data_dir + %d%central_atom_id
	
	path_to_input_files = "/Users/ltian/ART_sample/Cu64Zr36/10E11/"
	path_to_root_dir = path_to_input_files + 'art'
	total_energy = -46021.0583322591
	box_dim = [54.2654] * 3
	sample_name = "dump_10E11"
	set_up_input_files(path_to_data_dir, input_param)
	"""
	path_to_input_files = input_param['path_to_input_files']
	total_energy = input_param['total_energy']
	box_dim = input_param['box_dim']
	box_range = input_param['box_range']
	sample_name = input_param['sample_name']
	sample_type = input_param['sample_type']
	
	central_atom_list = input_param['central_atom_list']
	
	if not os.path.isdir(path_to_data_dir):
		os.makedirs(path_to_data_dir)	
	
	# obtain the refconfig file from the dump or lammps data file, need total energy, box_dim
	if sample_type == 'dump':
		dump_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim, box_range)
	elif sample_type == 'lammps_data':
		data_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim)
	else:
		raise Exception('the sample type can be either dump or lammps_data, please select the correct sample!')
	potential_file = raw_input("Please type the potential file name:")
	# manual: modify the in.lammps file 1) modify the dump file name of read_dump
	# 2) modify potential file string
	
	for central_atom in central_atom_list:
		path_to_run_dir = os.path.join(path_to_data_dir, str(central_atom))
		# create a directory for each central_atom id under path_to_data_dir = path_to_root_dir + sample_name
		if not os.path.isdir(path_to_run_dir):
			os.makedirs(path_to_run_dir)
		# modify the central_atom in bart.sh, with iteration from central_atom_list
		modify_bart(path_to_input_files, central_atom)
		# copy all modified files into the path_to_run_dir, mod_bart.sh,
		# src_files = os.listdir(path_to_input_files)
		#for file_name in src_files:
		for file_name in ['mod_bart.sh','refconfig',input_param['sample_name'],"in.lammps",potential_file]:
		#if file_name == 'bart.sh' or file_name == "central_atom_list.json" or file_name == "interested_atom_list.json":
		#	continue
			full_file_name = os.path.join(path_to_input_files, file_name)
			if (os.path.isfile(full_file_name)):
				shutil.copy(full_file_name, path_to_run_dir)
	print "finished setting up ART run dir and input files!"
	
	
def modify_bart(path_to_input_files, central_atom_id):
	"""
	this function takes the bart.sh in path_to_input_files and modify the central_atom
	then return the modified string file to be ready saved into a file 
	"""
	path_to_bart = os.path.join(path_to_input_files,'bart.sh')
	with open(path_to_bart, 'r') as f:
		data = f.read()
	pattern = "(setenv[\s]+Central_Atom[\s]+)([\d]+)"
	modified_bart = re.sub(pattern, r"\g<1>%s"%str(central_atom_id), data)
	
	path_to_mod_bart = os.path.join(path_to_input_files,'mod_bart.sh')
	with open(path_to_mod_bart, 'w') as f:
		f.write(modified_bart)
		f.close()
    
def data_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim):
	"""
	this function takes the data in the lammps data file and converts into 
	a refconfig file with periodic boundary condition
	refconfig file examples all used the absolute coordinates
	This is preferred way
	"""
	path_to_data = os.path.join(path_to_input_files, sample_name)
	
	result = read_data_from_lammps_data(path_to_data)
	
	df = result[['atom_id', 'x','y','z']]
	df.insert(loc=0, column='1', value='')
	file_to_save = os.path.join(path_to_input_files,'refconfig')
	with open(file_to_save, 'w') as f:
		f.write(' run_id: 1000\n total_energy: %s \n P %s %s %s\n'%(total_energy, box_dim[0],box_dim[1],box_dim[2]))
		df.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving refconfig in %s"%path_to_input_files	

def dump_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim, box_range, fractional_coordinates=True):
	"""
	this function takes a dump file and modify the dump file to a refconfig
	file for periodic boundary condition
	path_to_input_files = "/Users/ltian/ART_sample/Cu64Zr36/10E11/"
	sample_name= "dump_10E11"
	total_energy = -46021.0583322591
	box_dim = [54.2654] * 3
	dump_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim)
	"""
	path_to_dump = os.path.join(path_to_input_files,sample_name)
	result = read_data_from_dump(path_to_dump)
	df = result[['atom_id', 'x','y','z']]
	df.insert(loc=0, column='1', value='')
	# the atomic coordinates in dump file are fractional coordinates
	# the actual coordinates are (x_s * box_dim[0]) + box_range_x_start,...
	if fractional_coordinates == True:
		df.update((df[['x','y','z']] * box_dim) + [box_range[0][0],box_range[1][0],box_range[2][0]])
	file_to_save = os.path.join(path_to_input_files,'refconfig')
	with open(file_to_save, 'w') as f:
		f.write(' run_id: 1000\n total_energy: %s \n P %s %s %s\n'%(total_energy, box_dim[0],box_dim[1],box_dim[2]))
		df.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving refconfig in %s"%path_to_input_files

def check_tests_status(path_to_data_dir, input_param):
	"""
	check the status of each test of central_atom_list in input_param under path_to_data_dir
	as finished test or unfinished test
	
	This function will always overwrite the most update test status result
	into input SETTINGS file input_tests_done.json and input_tests_undone.json
	for the following processing
	
	Return:
		tests_undone: list
			the list of test id that has not been done yet
	"""
	
	# read the Max_num_events from bart.sh in path_to_input_files
	path_to_input_files = input_param['path_to_input_files']
	path_to_bart = os.path.join(path_to_input_files,"bart.sh")
	
	with open(path_to_bart, 'r') as f:
		data = f.read()
	pattern = "(setenv[\s]+Max_Number_Events[\s]+)([\d]+)"
	match = re.search(pattern, data)
	max_num_events = int(match.group(2))
	
	#list_of_test_id = input_param["list_of_test_id"]
	central_atom_list = input_param['central_atom_list']
	tests_not_done = []
	tests_done = []
	for test in central_atom_list:
		try:
			path_to_test = data_dir_to_test_dir(path_to_data_dir, test)
		except Exception:
			tests_not_done.append(test)
		# check if this test contains the final configuration file
		# based on max_num_events
		final_min_id = 1000 + max_num_events
		path_to_final_min = os.path.join(path_to_test,"min%s"%final_min_id)
		if os.path.isfile(path_to_final_min):
			tests_done.append(test)
		else:
			tests_not_done.append(test)
	print "In %s , finished tests ids in central_atom_list of current input SETTINGS file are:"%path_to_data_dir, tests_done
	print "In %s , un-finished tests ids in central_atom_list of current input SETTINGS file are:"%path_to_data_dir, tests_not_done
	print "In %s, start creating art_data input SETTINGs files for finished tests (input_tests_done.json) and unfinished tests(input_tests_undone.json):"%path_to_data_dir
	
	path_to_tests_done = os.path.join(path_to_data_dir,"input_tests_done.json")
	path_to_tests_undone = os.path.join(path_to_data_dir,"input_tests_undone.json")
	
	input_param_tests_done = copy.deepcopy(input_param)
	input_param_tests_undone = copy.deepcopy(input_param)
	input_param_tests_done["central_atom_list"] = tests_done
	input_param_tests_undone["central_atom_list"] = tests_not_done
	
	input_param_tests_done["list_of_test_id"] = tests_done
	input_param_tests_undone["list_of_test_id"] = tests_not_done
	
	with open(path_to_tests_done, 'w') as f:
		json.dump(input_param_tests_done,f,indent=2)
	with open(path_to_tests_undone, 'w') as f:
		json.dump(input_param_tests_undone,f,indent=2)
	
	print "\n"
	print "For finished tests:"
	print "Now user can check only the results of finished tests by using art_data -s input_tests_done.json --filter, art_data -s input_tests_done.json --eng --calc etc"
	print "For unfinished tests:"
	print "Now user can choose to delete these unfinished tests completely by using art_data -s input_tests_undone.json --art --delete_tests"
	print "Then user can continue to run these unfinished tests from the beginning by using art_data -s input_tests_undone.json --art --run"
	print "test status check done!"
	return tests_not_done
	
def delete_art_tests_files(path_to_data_dir, file_names):
	"""
	file_names is a list of files to be deleted in the final art tests directory
	"""
	path_to_central_atom_list = os.path.join(path_to_data_dir,"central_atom_list.json")
	central_atom_list = json.load(open(path_to_central_atom_list, 'r'))
	i=0
	for central_atom in central_atom_list:
		for each_file in file_names:
			path_to_file = os.path.join(path_to_data_dir, str(central_atom), each_file)
			if os.path.isfile(path_to_file):
				if i >= 1:
					print "deleting the file %s"%path_to_file
					os.remove(path_to_file)
				else:
					print "deleting the file %s"%path_to_file
					print "confirm deleting (y/n):"
					if prompt_yes_no() is True:
						os.remove(path_to_file)
			else:
				print "%s"%each_file, "does not exist in %s"%os.path.join(path_to_data_dir, str(central_atom))
		i = i + 1

def delete_art_tests(path_to_data_dir, central_atom_list):
	"""
	file_names is a list of files to be deleted in the final art tests directory
	"""
	print "begining delete tests:", central_atom_list
	print "confirm deleting (y/n):"
	if prompt_yes_no() is True:
		for test in central_atom_list:
			try:
				path_to_test = data_dir_to_test_dir(path_to_data_dir, test)
				if os.path.isdir(path_to_test):
					print "deleting the test %s"%path_to_test
					shutil.rmtree(path_to_test)
				else:
					print "test %s"%test, "does not exist in %s"%path_to_data_dir
			except Exception:
				print "test %s"%test, "does not exist in %s"%path_to_data_dir
	else:
		print "Not deleting!"
	print "done deleting all tests in specified list of tests in %s"%path_to_data_dir

def delete_unused_events_data(path_to_data_dir,input_param):
	"""
	this function delete the configuration files of un-used events data
	inside each test of an art data project to save disk space, especially
	when users has not used WRITE_REJECTED_EVENTS = .False. option in bart.sh
	
	First, it will read the final_selected_events.json, which stores all useful events.
	
	Second, it will delete all configuration files that are not saved
	in the final_selected_events.json
	"""
	path_to_final_selected_events = os.path.join(path_to_data_dir,"final_selected_events.json")
	if os.path.isfile(path_to_final_selected_events):
		print "reading final_selected_events.json, ensure that you always get most updated final_selected_events.json by --filter --re_calc if you have calculated more tests"
		final_selected_events = json.load(open(path_to_final_selected_events,"r"))
	else:
		raise Exception("final_selected_events.json does not exist in %s"%path_to_data_dir)

	all_tests_events = dict()
	for event in final_selected_events:
		test_id = int(event[0][4:])
		init, sad, fin = event[1][0],event[1][1],event[1][2]
		if test_id in all_tests_events:
			all_tests_events[test_id].extend([init,sad,fin])
		else:
			all_tests_events[test_id] = [init,sad,fin]
			#all_tests_id.append(test_id)
	
	print "confirm deleting (y/n):"
	if prompt_yes_no() is False:
		print "response received, not deleting"
		return None
	else:
		print "response received, begin deleting"
	
	path_to_central_atom_list = os.path.join(path_to_data_dir,"central_atom_list.json")
	if os.path.isfile(path_to_central_atom_list) or "central_atom_list" in input_param:
		if os.path.isfile(path_to_central_atom_list):
			print "reading central_atom_list.json"
			central_atom_list = json.load(open(path_to_central_atom_list, 'r'))	
		elif "central_atom_list" in input_param:
			print "reading central_atom_list from input SETTINGs file"
			central_atom_list = input_param['central_atom_list']
		delete_tests_list = []
		saved_tests_list = []
		for test_id in central_atom_list:
			if test_id not in all_tests_events:
				delete_tests_list.append(test_id)
			else:
				saved_tests_list.append(test_id)
		delete_art_tests(path_to_data_dir,delete_tests_list)
		print ">>> confirm updating central_atom_list.json file (y/n): save original central_atom_list.json if necessary"
		if prompt_yes_no() is True:
			print "response received, updating!"
			update_central_atom_list(path_to_data_dir,saved_tests_list)
		else:
			print "response received, not updating"
	else:
		print "central_atom_list.json does not exist in %s"%path_to_data_dir
		print "central_atom_list key does not exist in input SETTING file"
		print ">>> only deleting unused events for all tests stored in final_selected_events.json"
	
	print "\n Now begin deleting unused events configuration files for all tests stored in final_selected_events.json"
	for test_id in all_tests_events:
		path_to_test_dir = data_dir_to_test_dir(path_to_data_dir,test_id)
		for f in os.listdir(path_to_test_dir):
			is_match_config = re.match(r"min[0-9]+", f) or re.match(r"sad[0-9]+", f) or re.match(r"min[0-9]+.dump", f) or re.match(r"sad[0-9]+.dump", f)
			path_to_file = os.path.join(path_to_test_dir, f)
			is_file = os.path.isfile(path_to_file)
			if is_match_config and is_file:
				if f.endswith('.dump'):
					config_id = f[:-5]
				else:
					config_id = f
				if config_id not in all_tests_events[test_id]:
					print "deleting the file %s"%path_to_file
					os.remove(path_to_file)
	print "done deleting unused events data!"

