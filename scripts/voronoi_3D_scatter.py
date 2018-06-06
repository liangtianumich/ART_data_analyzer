#!/usr/bin/env python
import os
from visualizer.voronoi_visualizer import voronoi_scatter_3d

path_to_data_dir = os.environ['DATA_DIR']

path_to_curr_event = path_to_data_dir + "test1/results/event_min1000_sad1001_min1001"

path_to_config = path_to_data_dir + "test1/" + "min1000.dump"

voronoi_scatter_3d(path_to_curr_event, path_to_config)
