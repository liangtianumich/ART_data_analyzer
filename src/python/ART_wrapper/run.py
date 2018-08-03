"""
this module runs ART simulations by set up the ART input files for each MD dump file sample, 
create a run_dir associated with the central_atom id, 
copy all modified input files into each run_dir, run ./bart.sh under each run_dir
by subprocess module in parallel
"""
import subprocess
import data_reader

def run_art(path_to_root_dir, sample_name, path_to_input_files, pbs=False):
	"""
	
	"""

def set_up_input_files(path_to_root_dir, sample_name, path_to_input_files):
	"""
	this function will set up the input files for each MD dump file sample,
	, create a run_dir associated with the central_atom id,
	copy all modified input files into each run_dir,
	run bart.sh under each run_dir
	
	directory structure: 
	
	# path_to_input_files are saved in the $ART_INPUT
	# path_to_data_dir = path_to_root_dir + sample_name
	# path_to_run_dir = path_to_data_dir + %d%central_atom_id
	"""
	
	
	# get the dump file that has string called sample name
	
	
	# obtain the refconfig file from the dump file, need total energy
	
	
	# modify the in.lammps file 1) modify the dump file name of read_dump
	# 2) modify potential file string
	
	
	# modify the central_atom in bart.sh, with iteration from sample dump file atom_id
	
	
	# create a directory for each central_atom id under path_to_data_dir = path_to_root_dir + sample_name
	
	# copy all modified files into the path_to_run_dir = path_to_data_dir + %d%central_atom_id
	
	# run bart.sh inside each dir
	



def dump_to_refconfig(path_to_input_files, sample_name):
	"""
	this function takes a dump file and modify the dump file to a refconfig
	file
	"""
	path_to_dump = path_to_input_files + sample_name
	result = read_data_from_dump(path_to_dump)
	
	
	
	
	
	


	
	
	
