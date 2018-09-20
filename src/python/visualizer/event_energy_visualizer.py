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

def plot_act_relax_histogram(path_to_image, x):
	"""
	path to image now point to the event dir
	x is a two element list, each element is an array
	"""
	plt.figure()
	weights = [np.ones_like(x[0])/float(len(x[0])), np.ones_like(x[1])/float(len(x[1]))]
	plt.hist(x,bins=10, weights=weights, color=['r','black'],label=["activation energy",'relaxation energy'])
	plt.xlabel("activation or relaxation energy /eV",fontsize=20)
	plt.ylabel("percentage", fontsize=20)
	plt.legend(loc='best')
	plt.savefig(path_to_image)
	plt.close()
