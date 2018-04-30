#from "" import *
#from "../src/python/data_reader" import *
import imp
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
strain_calculator = imp.load_source('module.name', '../src/python/calculator/strain_calculator.py')
data_reader = imp.load_source('module.name', '../src/python/data_reader.py')
strain_visualizer = imp.load_source('module.name', '../src/python/visualizer/strain_visualizer.py')

path_to_file = "strain_result.pkl"
path_to_ovito_file_xyz = "ovito_strain_xyz.txt"

if os.path.exists(path_to_file):
	strain_results = pickle.load(open(path_to_file,"rb"))
	ovito_strain_xyz = pd.read_csv(path_to_ovito_file_xyz, sep='\s+', skiprows = 1)
	ovito_strain_xyz = ovito_strain_xyz.iloc[:,5:7]
	#ovito_strain.columns = ["vol_strain","shear_strain","item_id","atom_type","x","y","z"]
	ovito_strain_xyz.columns = ["vol_strain","shear_strain"]
	
	atom_item = []
	shear_strain = []
	ovito_xyz = []
	for index,row in ovito_strain_xyz.iterrows():
		ovito_xyz.append(row["shear_strain"])
		atom_item.append(index + 1)
		#row_xyz = ovito_strain_xyz.iloc[index]
		#ovito_xyz.append(row_xyz["shear_strain"])
		shear_strain.append(strain_results[index+1][1])
	
	plt.plot(atom_item,ovito_xyz,'ro',label="ovito_xyz")
	plt.plot(atom_item,shear_strain,'bx',label="my code")
	plt.xlabel("atom item id")
	plt.ylabel("shear strain")
	plt.legend(loc='best')
	plt.savefig("shear_strain_comparison.tif")
	
else:
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


