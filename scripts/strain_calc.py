import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from data_reader import *
from calculator.strain_calculator import *
from visualizer.strain_visualizer import *
from util import Atom


path_to_test_dir = os.environ['TEST_DIR']
path_to_input = path_to_test_dir + "/input.json"
path_to_results = path_to_test_dir + "/results.pkl"

if os.path.exists(path_to_results):
	strain_results = pickle.load(open(path_to_results,"rb"))
	atom_item = []
	shear_strain = []
	for (key,pair) in strain_results.items():
		atom_item.append(key)
		shear_strain.append(pair[1])
	plt.plot(atom_item,shear_strain,'bx',label="shear_strain")
	plt.xlabel("atom item id")
	plt.ylabel("shear strain")
	plt.legend(loc='best')
	plt.savefig("shear_strain.tif")
	
else:
	path_to_example = os.environ['MY_ART'] + "/examples/example_1/"
	
	path_to_file_ini = path_to_example + "min1000.dump"
	path_to_file_sad = path_to_example + "min1001.dump"
	
	initial_config_data = read_data_from_dump(path_to_file_ini)
	saddle_config_data = read_data_from_dump(path_to_file_sad)

	cut_off_distance = {(1,1):3.7,(1,2):3.7,(2,2):3.7}

	size = 32.130125 - 0.299875
	box_dim = [size, size, size]
	
	calc_input = {'cut_off': cut_off_distance, 'box_dim': box_dim}
	if not os.path.exists(path_to_test_dir):
		os.makedirs(path_to_test_dir)
	
	with open(path_to_input,'w') as f:
		pickle.dump(calc_input, f)
		f.close()
		
	
	strain_result = local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, path_to_results)
	
	
	x, y, z = [],[],[]
	color_value = []
	for key, atom_strain in strain_result.items():	
		atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
		x.append((atom.atom_loc)[0])
		y.append((atom.atom_loc)[1])
		z.append((atom.atom_loc)[2])
		color_value.append(atom_strain[0])
	scatter3d(x,y,z, color_value, colorsMap='jet')


