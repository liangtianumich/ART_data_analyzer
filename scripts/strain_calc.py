#from "" import *
#from "../src/python/data_reader" import *
import imp
strain_calculator = imp.load_source('module.name', '../src/python/calculator/strain_calculator.py')
data_reader = imp.load_source('module.name', '../src/python/data_reader.py')
strain_visualizer = imp.load_source('module.name', '../src/python/visualizer/strain_visualizer.py')


path_to_file_ini = "/Users/ltian/Documents/alabama/ART_data_analyzer/examples/example_1/min1000.dump"
initial_config_data = data_reader.read_data_from_dump(path_to_file_ini)

path_to_file_sad = "/Users/ltian/Documents/alabama/ART_data_analyzer/examples/example_1/min1001.dump"
saddle_config_data = data_reader.read_data_from_dump(path_to_file_sad)

cut_off_distance = {(1,1):3.7,(1,2):3.7,(2,2):3.7}

size = 32.130125 - 0.299875
box_dim = [size, size, size]

strain_result = strain_calculator.local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim)
x, y, z = [],[],[]
color_value = []
for key, atom_strain in strain_result.items():	
	atom = strain_calculator.Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
	x.append((atom.atom_loc)[0])
	y.append((atom.atom_loc)[1])
	z.append((atom.atom_loc)[2])
	color_value.append(atom_strain[0])
strain_visualizer.scatter3d(x,y,z, color_value, colorsMap='jet')

