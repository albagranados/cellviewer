# Date: 01 / 2018

# modules
import sys
import os
import numpy as np
import matplotlib

# my_modules
sys.path.append(os.path.join(os.getcwd(), "src"))  # if we want to run in the Python Console
import iprocessing as iproc  # image processings
import utilities as util

# configure latex plots
matplotlib.rcParams['text.usetex'] = True; matplotlib.rcParams['font.family'] = 'serif'

# module info
__author__ = 'alba granados'
__email__ = "alba.granados@crg.eu"
__project__ = 'voronoi analysis'

# set working directory
parent_dir = '/home/alba/DropboxCRG/postdoc_CRG/coding/cellviewer/voronoi'
os.chdir(parent_dir)
# set output directory
output_dir = parent_dir + '/output/'
if not os.path.exists(output_dir):  os.makedirs(output_dir)

# # ========= GENERATE SYNTHETIC POINT PATTERN
# # ==================
dict_synthetic_pp = {'unit_size': 1, 'cell_size': 20 ** 2, 'n_points_total': 1000,
                       'n_clusters': 1, 'distribution': 'uniform',
                       'geometry': 'circle', 'enrichment_ratio': 40,
                       'mean_radius': 0.5, 'std_radius': 0,
                       'radius': np.append(np.arange(20, 80, 5), np.matlib.repmat(50, 1, 10))}
#    n_cb = 0.5
#    dict_synthetic_data['enrichment_ratio'] = n_cb*((dict_synthetic_data.get('cell_size')-dict_synthetic_data.get(
#        'n_clusters')*np.pi*dict_synthetic_data.get('mean_radius')**2)/(np.pi*dict_synthetic_data.get(
#        'mean_radius')**2))
#    n_points_clusters = 60
#    dict_synthetic_data['n_points_total'] = n_points_clusters*(dict_synthetic_data.get('n_clusters')+1./dict_synthetic_data.get('enrichment_ratio')*((dict_synthetic_data.get('cell_size')-dict_synthetic_data.get(
#        'n_clusters')*np.pi*dict_synthetic_data.get('mean_radius')**2)/(np.pi*dict_synthetic_data.get(
#        'mean_radius')**2)))
points = util.generate_pattern_rnd(dict_synthetic_pp, num=1, save_data={'save': 1, 'output_dir': output_dir})

# # ============ CROP & SAVE
# # ==================
crop = [168, 176, 75, 79]
points11 = iproc.points_2dcrop(points[np.where(channels == 1)], crop[0], crop[1], crop[2], crop[3])
np.savetxt(output_dir + 'test1.txt', points11)
points22 = iproc.points_2dcrop(points[np.where(channels == 2)], crop[0], crop[1], crop[2], crop[3])
np.savetxt(output_dir + 'test2.txt', points22)