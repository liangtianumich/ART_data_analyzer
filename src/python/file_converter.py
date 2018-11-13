import os
from data_reader import read_data_from_non_dump
from util import data_dir_to_test_dir

def refconfig_to_lammps_data(path_to_data_dir, event_state, input_param):
	"""
	this function converts the configuration files of a single event stored
	in event_state from ART output refconfig type into lammps data type
	"""	
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	init, sad, fin = event_state[1][0],event_state[1][1],event_state[1][2]
	#path_to_event = path_to_test_dir + "/results/event_" + init + "_" + sad + "_" + fin
	
	path_to_file_ini = path_to_test_dir + '/' + init
	path_to_file_sad = path_to_test_dir + '/' + sad
	path_to_file_fin = path_to_test_dir + '/' + fin
	
	ini_config_data = read_data_from_non_dump(path_to_file_ini)
	sad_config_data = read_data_from_non_dump(path_to_file_sad)
	fin_config_data = read_data_from_non_dump(path_to_file_fin)
	
	ini_config_data = ini_config_data[["item","atom_id","x","y","z"]]
	sad_config_data = sad_config_data[["item","atom_id","x","y","z"]]
	fin_config_data = fin_config_data[["item","atom_id","x","y","z"]]
	
	num_of_atoms = ini_config_data.shape[0]
	item_list = [x+1 for x in range(num_of_atoms)]
	
	num_atom_types = len(ini_config_data['atom_id'].unique())
	#print "ini_config_data:", ini_config_data
	#ini_config_data.insert(loc=0, column='item', value=item_list)
	#sad_config_data.insert(loc=0, column='item', value=item_list)
	#fin_config_data.insert(loc=0, column='item', value=item_list)
	
	ini_file_to_save = os.path.join(path_to_test_dir, init+".lammps")
	sad_file_to_save = os.path.join(path_to_test_dir, sad+".lammps")
	fin_file_to_save = os.path.join(path_to_test_dir, fin+".lammps")
	
	box_range = input_param["box_range"]
	
	with open(ini_file_to_save, 'w') as f:
		f.write('# LAMMPS data file written by OVITO\n%s atoms\n%s atom types\n%s %s xlo xhi\n%s %s ylo yhi\n%s %s zlo zhi\n\nAtoms # atomic\n\n'%(num_of_atoms,num_atom_types,box_range[0][0],box_range[0][1],box_range[1][0],box_range[1][1],box_range[2][0],box_range[2][1]))
		ini_config_data.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving file in %s"%ini_file_to_save
	
	with open(sad_file_to_save, 'w') as f:
		f.write('# LAMMPS data file written by OVITO\n%s atoms\n%s atom types\n%s %s xlo xhi\n%s %s ylo yhi\n%s %s zlo zhi\n\nAtoms # atomic\n\n'%(num_of_atoms,num_atom_types,box_range[0][0],box_range[0][1],box_range[1][0],box_range[1][1],box_range[2][0],box_range[2][1]))
		sad_config_data.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving file in %s"%sad_file_to_save
	
	with open(fin_file_to_save, 'w') as f:
		f.write('# LAMMPS data file written by OVITO\n%s atoms\n%s atom types\n%s %s xlo xhi\n%s %s ylo yhi\n%s %s zlo zhi\n\nAtoms # atomic\n\n'%(num_of_atoms,num_atom_types,box_range[0][0],box_range[0][1],box_range[1][0],box_range[1][1],box_range[2][0],box_range[2][1]))
		fin_config_data.to_csv(f, header=False, index=False, sep=' ')
		f.close()
	print "saving file in %s"%fin_file_to_save
	
	
	
	
