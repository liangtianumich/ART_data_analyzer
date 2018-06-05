"""
this vornoi visualization module plot the statistical distribution of voronoi classes
i.e, ico, ico-like, GUMs
"""
import os
import json
import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

def plot_voronoi_histogram_3(path_to_image, x):
	"""
	Input:
	path_to_image: string
		string that specifies the path to save the image
	x: list
		a three element list, each element is an array-like
	
	"""
	plt.figure()
	#fig, ax = plt.subplots()
	#weights = [np.ones_like(x[0])/float(len(x[0])), np.ones_like(x[1])/float(len(x[1])), np.ones_like(x[2])/float(len(x[2]))]
	plt.hist(x, bins=10, color=['r','b','black'],label=["initial","saddle","final"])
	x1 = [0,1,2]
	labels = ["ico","ico-like","GUM"]
	#ax.set_xticks(x1)
	#ax.set_xticklabels(labels)
	plt.xticks(x1, labels)
	plt.legend(loc='best')
	plt.savefig(path_to_image)
	plt.close()
