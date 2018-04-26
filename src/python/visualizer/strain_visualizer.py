"""
this strain visualization module plot the distribution of local atomic strains
in terms of atomic coordinates x,y,z in the simulation box
"""
import pandas as pd
from strain_calculator import Atom
import matlibplot.pyplot as plt

def plot_strain_2d(initial_config_data, strain, normal_plane):
	"""
	
	this function plot the distributions of local atomic strains
	in 2d for von mises strain and hydrostatic invariants
	
	"""
	for key,atom_strain in strain.items():
		
		atom = Atom.from_ds(initial_config_data.loc[initial_config_data["item"] == key])
		atom_strain[0]
		atom_strain[1]
		
