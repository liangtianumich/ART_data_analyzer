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
from data_reader import *
import shutil


def run_art(path_to_root_dir, sample_name, path_to_input_files, pbs=False):
	"""
	
	"""
	# run bart.sh inside each dir

def set_up_input_files(path_to_root_dir, path_to_input_files, total_energy, box_dim, sample_name, sample_type='dump', path_to_central_atom_list=None):
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
	set_up_input_files(path_to_root_dir, path_to_input_files, total_energy, box_dim, sample_name, sample_type='dump', path_to_central_atom_list=None)
	"""
	
	
	# get the path string of dump or lammps data file
	path_to_sample = os.path.join(path_to_input_files,sample_name)
	path_to_data_dir = os.path.join(path_to_root_dir, sample_name)
	if not os.path.isdir(path_to_data_dir):
		os.makedirs(path_to_data_dir)
		
	
	# obtain the refconfig file from the dump or lammps data file, need total energy, box_dim
	if sample_type == 'dump':
		dump_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim)
	elif sample_type == 'lammps_data':
		data_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim)
	else:
		raise Exception('the sample type can be either dump or lammps_data, please select the correct sample!')
	
	# manual: modify the in.lammps file 1) modify the dump file name of read_dump
	# 2) modify potential file string
	
	
	# modify the central_atom in bart.sh, with iteration from sample dump file atom item id
	if path_to_central_atom_list is None:
		
		path_to_central_atom_list = os.path.join(path_to_data_dir, "central_atom_list.json")
		
		if sample_type == 'dump':
			config_results = read_data_from_dump(path_to_sample)
		elif sample_type == 'lammps_data':
			config_results = read_data_from_lammps_data(path_to_sample)
		list_of_test_id = config_results['item'].tolist()
		
		with open(path_to_central_atom_list, 'w+') as f:
			json.dump(list_of_test_id, f)
	elif os.path.isfile(path_to_central_atom_list):
		list_of_test_id = json.load(open(path_to_central_atom_list, 'r'))
	else:
		raise Exception("%s is not a file"%path_to_central_atom_list)
	
	for central_atom in list_of_test_id:
		path_to_run_dir = os.path.join(path_to_data_dir, str(central_atom))
		# create a directory for each central_atom id under path_to_data_dir = path_to_root_dir + sample_name
		if not os.path.isdir(path_to_run_dir):
			os.makedirs(path_to_run_dir)
		
		modify_bart(path_to_input_files, central_atom)
		# copy all modified files into the path_to_run_dir, mod_bart.sh,
		src_files = os.listdir(path_to_input_files)
		for file_name in src_files:
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

def dump_to_refconfig(path_to_input_files, sample_name, total_energy, box_dim):
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
	file_to_save = os.path.join(path_to_input_files,'refconfig')
	with open(file_to_save, 'w') as f:
		f.write(' run_id: 1000\n total_energy: %s \n P %s %s %s\n'%(total_energy, box_dim[0],box_dim[1],box_dim[2]))
		df.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving refconfig in %s"%path_to_input_files
