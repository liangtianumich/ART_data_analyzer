import os
import json
import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from util import Atom

def plot_2d(path_to_dir, x, y, x_label, y_label):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(x, y,'ro',markersize=1.5)
	plt.xlabel(x_label,fontsize=20)
	plt.ylabel(y_label,fontsize=20)
	plt.savefig(path_to_dir)
	plt.close()

def plot_2d_train_fit(path_to_dir, x, y,x_fit,y_fit, x_label, y_label):
	
	"""
	this function plot a 2d plot for displacement and strain for a single event
	
	displacement:list
	
	strain:list
	"""
	plt.figure()
	plt.plot(x, y,'rX',markersize=2.0)
	plt.plot(x_fit, y_fit,'b')
	plt.xlabel(x_label,fontsize=18)
	plt.ylabel(y_label,fontsize=18)
	plt.savefig(path_to_dir)
	plt.close()

def plot_histogram_3(path_to_image, x):
	"""
	path to image now point to the event dir
	x is a three element list, each element is an array
	"""
	
	plt.figure()
	#remove nan elements
	new_x = []
	for item in x:
		item_array = np.array(item)
		new_x.append(item_array[~np.isnan(item_array)].tolist())
	x=new_x
	
	weights = [np.ones_like(x[0])/float(len(x[0])), np.ones_like(x[1])/float(len(x[1])), np.ones_like(x[2])/float(len(x[2]))]
	plt.hist(x, bins=10, weights=weights, color=['r','b','black'], label=["initial-saddle","saddle-final","initial-final"], edgecolor='black', linewidth=1.2)
	plt.legend(loc='best')
	plt.savefig(path_to_image)
	plt.close()

def plot_histogram(path_to_image, x, color='b'):
	plt.figure()
	weights = np.ones_like(x)/float(len(x))
	plt.hist(x, bins=30, weights=weights, color=color, edgecolor='black', linewidth=1.2)
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
	scatter3d(x,y,z, color_value)

def scatter3d(x,y,z, cs, str_name="strain", colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig('%s_3d.png'%str_name)
    plt.close()

def scatter3d_path(path_to_image, x,y,z,cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    # view the 3D plot in x-y plane
    ax.view_init(azim=0, elev=90)
    # 2D projection whose z =0 with zdir is z
    ax.scatter(x, y, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig(path_to_image,dpi=600)
    plt.close()
