"""
this strain visualization module plot the distribution of local atomic strains
in terms of atomic coordinates x,y,z in the simulation box
"""
import pandas as pd
from util import Atom
import matplotlib
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
		
