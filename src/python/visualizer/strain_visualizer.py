"""
this strain visualization module plot the distribution of local atomic strains
in terms of atomic coordinates x,y,z in the simulation box
"""
import os
import pickle
import pandas as pd
from util import Atom
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
# install mpl_toolkits module may need to add the __init__.py manually in
# site-packages/ to make it a package, if installed mpl_toolkits
# if mpl_tookits package not included along with matplotlib, sudo pip install -U matplotlib
# should do the trick, if not, i tried brew install geos, then sudo pip install - U matplotlib,
# it upgrades successfully, then sys.path updated and importlib.import_module('mpl_toolkits').__path__
# points to the symlink, brew remove python@2 --ignore-dependencies
# to remove the python2 dependencies created by brew to use system default python dependencies
def strain_events_stats_visualization(path_to_data_dir, list_of_test_id):
	"""
	this function visualize the strain statistics from tests listed in
	list_of_test_id, the user can customize the tests they want to know
	their statistics, among these tests, only the one has events_stat.pkl
	file will be used for plot
	
	path_to_data_dir: str
		path to the data dir
	list_of_tests: list
		a number list showing the test id number
		e.g. list [1,2] mean that visualizing test1, test2 data
		
	"""
	# set up the statistical quantities for the tests
	disp_ave, disp_std, disp_max , disp_ave_2, disp_std_2, disp_max_2 = [], [], [], [], [], []
	
	shear_ave, shear_std, shear_max, shear_ave_2, shear_std_2, shear_max_2 = [], [], [], [], [], []
	
	vol_ave, vol_std, vol_max, vol_ave_2, vol_std_2, vol_max_2 = [], [], [], [], [], []
	
	for i in list_of_test_id:
		path_to_curr_test = path_to_data_dir + "test%s"%i
		path_to_curr_event = path_to_curr_test + "/results/events_stats.pkl"
		# skip the test who do not have events_stats.pkl in all tests specified in
		# list_of_test_id
		if os.path.exists(path_to_curr_event):
			result = pickle.load(open(path_to_curr_event,'r'))
			for event in result:
				init_sad = event[0]
				sad_fin = event[1]
				# calculate the statistics of init_sad and sad_fin		
				disp_ave.append(init_sad["ave"][2])
				disp_std.append(init_sad["std"][2])
				disp_max.append(init_sad["max"][2])
				
				shear_ave.append(init_sad["ave"][1])
				shear_std.append(init_sad["std"][1])
				shear_max.append(init_sad["max"][1])
				
				vol_ave.append(init_sad["ave"][0])
				vol_std.append(init_sad["std"][0])
				vol_max.append(init_sad["max"][0])
				
				disp_ave_2.append(sad_fin["ave"][2])
				disp_std_2.append(sad_fin["std"][2])
				disp_max_2.append(sad_fin["max"][2])
				
				shear_ave_2.append(sad_fin["ave"][1])
				shear_std_2.append(sad_fin["std"][1])
				shear_max_2.append(sad_fin["max"][1])
				
				vol_ave_2.append(sad_fin["ave"][0])
				vol_std_2.append(sad_fin["std"][0])
				vol_max_2.append(sad_fin["max"][0])
	
	plot_histogram_2(path_to_data_dir+"/disp_ave.png", [disp_ave,disp_ave_2])
	plot_histogram_2(path_to_data_dir+"/disp_std.png", [disp_std,disp_std_2])
	plot_histogram_2(path_to_data_dir+"/disp_max.png", [disp_max,disp_max_2])
	
	plot_histogram_2(path_to_data_dir+"/shear_ave.png", [shear_ave,shear_ave_2])
	plot_histogram_2(path_to_data_dir+"/shear_std.png", [shear_std,shear_std_2])
	plot_histogram_2(path_to_data_dir+"/shear_max.png", [shear_max,shear_max_2])
	
	plot_histogram_2(path_to_data_dir+"/vol_ave.png", [vol_ave,vol_ave_2])
	plot_histogram_2(path_to_data_dir+"/vol_std.png", [vol_std,vol_std_2])
	plot_histogram_2(path_to_data_dir+"/vol_max.png", [vol_max,vol_max_2])	
	print "done plotting strain statistics for all interested tests!"

def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig('strain_3d.png')
    plt.close()

def plot_histogram(path_to_image, x):
	plt.figure()
	plt.hist(x,bins='auto')
	plt.savefig(path_to_image)
	plt.close()
	
def plot_histogram_2(path_to_image, x):
	"""
	path to image now point to the event dir
	x is a two element list, each element is an array
	"""
	
	plt.figure()
	weights = [np.ones_like(x[0])/float(len(x[0])), np.ones_like(x[1])/float(len(x[1]))]
	plt.hist(x,bins=10, weights=weights, color=['r','black'],label=["initial-saddle","saddle-final"])
	plt.legend(loc='best')
	plt.savefig(path_to_image)
	plt.close()
	
	
def scatter_plot_3D(strain_result, initial_config_data):
	"""
	this function plot atomic strain, shear strain in 3d
	strain_result is the output dict of strain_calculator.local_strain_calculator_orth
	it still need the initial configuration data to extract the atomic coordinates
	input arguments:
	strain_results: dict()
	
	initial_config_data: pandas.Dataframe
	"""
	x, y, z = [],[],[]
	color_value = []
	for key, atom_strain in strain_result.items():	
		atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
		x.append((atom.atom_loc)[0])
		y.append((atom.atom_loc)[1])
		z.append((atom.atom_loc)[2])
		color_value.append(atom_strain[0])
	scatter3d(x,y,z, color_value, colorsMap='jet')

def plot_2d_shear(path_to_dir, displacement, strain):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(displacement, strain,'ro')
	plt.xlabel('atomic displacement',fontsize=20)
	plt.ylabel('von Mises shear strain',fontsize=20)
	plt.savefig(path_to_dir)
	plt.close()

def plot_2d_vol(path_to_dir, displacement, strain):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(displacement, strain,'ro')
	plt.xlabel('atomic displacement',fontsize=20)
	plt.ylabel('volumetric strain',fontsize=20)
	plt.savefig(path_to_dir)
	plt.close()
	
def plot_strain_2d(initial_config_data, strain, normal_plane):
	"""
	
	this function plot the distributions of local atomic strains
	in 2d for von mises strain and hydrostatic invariants
	
	"""
	for key,atom_strain in strain.items():
		
		atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
		atom_strain[0]
		atom_strain[1]
		
