# modules
import sys, os, time
import numpy as np, matplotlib.pyplot as plt
import cv2
from scipy.spatial import Voronoi
# my_modules
from src import vprocessing as vproc, iprocessing as iproc, statistics as stat, utilities as util

# def main():
#     # all the stuff

reload(stat)
plt.close('all')
# set directories
parent_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis'   # os.chdir(parent_dir)
source_dir = parent_dir + '/src/'; sys.path.insert(0, source_dir)
output_dir = parent_dir + '/output/'  # set output directory
if not os.path.exists(output_dir): os.makedirs(output_dir)

# # =========== READ INPUT FILE =================
# # ================================================
experiment_author = ''; file_dir = ''; file_name = ''

# # ====== MULTIPLE IMAGE STATISTICS =========
# # ==========================================
file_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/test/synthetic_pp/validation/6/output/'
fileExt = '.txt'

num_clusters = []; cluster_diameters = []; image_area = []; num_loc_cluster = []; ntotal = []
for fileid in os.listdir(file_dir):
    file_name = fileid.split(".txt")[0]
    path = os.path.join(file_dir, file_name + fileExt)

    with open(path) as f:
        for line in f:
            if line.startswith('number_clusters'):
                num_clusters.append(int(line.split('\t')[1]))
            elif line.startswith('cluster_diameters'):
                aux = line.split('\t')[1]
                aux = aux.replace('[', '').replace(']', '').split(',')
                aux = [float(ii) for ii in aux]
                cluster_diameters.append(aux)
            elif line.startswith('image_area'):
                image_area.append(float(line.split('\t')[1]))
            elif line.startswith('ntotal'):
                ntotal.append(int(line.split('\t')[1]))
            elif line.startswith('number_localizations_cluster'):
                aux = line.split('\t')[1]
                aux = aux.replace('[', '').replace(']', '').split(',')
                aux = [float(ii) for ii in aux]
                num_loc_cluster.append(aux)

cluster_diameters = [item for blob in cluster_diameters for item in blob]
percentage_num_loc = 100*np.divide(np.asarray([np.sum(loc_in_image) for loc_in_image in
                                                num_loc_cluster]), np.asarray(ntotal))
num_loc_cluster = [item for blob in num_loc_cluster for item in blob]

file_name = "_".join(file_name.split('_')[0:-2]) + '_statistics'

reload(stat)
stat.plot_hist(np.asarray(num_clusters), num_bins=20, xlabel=r'number of clusters [nm]')
plt.savefig(output_dir + file_name + '_numclusters_hist.pdf', bbox_inches='tight')

stat.plot_hist(0.5*np.asarray(cluster_diameters), num_bins=30, xlabel=r'cluster radius R [nm]')
plt.savefig(output_dir + file_name + '_clusterradius_hist.pdf', bbox_inches='tight')

stat.plot_hist(np.asarray(num_loc_cluster), num_bins=20, xlabel=r'N$^{cluster}$ [points/cluster]')
plt.savefig(output_dir + file_name + '_clusternumloc_hist.pdf', bbox_inches='tight')

stat.plot_hist(percentage_num_loc, num_bins=20, xlabel=r'percentage of localizations in clusters$ [$\%$]')
plt.savefig(output_dir + file_name + '_percloccluster_hist.pdf', bbox_inches='tight')

stat.plot_hist(np.divide(np.asarray(num_clusters), np.asarray(image_area)),
               num_bins=20, xlabel=r'$\kappa$ [clusters/nm$^2$]')
plt.savefig(output_dir + file_name + '_densityclusters_hist.pdf', bbox_inches='tight')

stat.plot_hist(np.divide(np.asarray(num_loc_cluster), (0.5*np.asarray(cluster_diameters))**2*np.pi),
               num_bins=20, xlabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
plt.savefig(output_dir + file_name + '_clusterdensities_hist.pdf', bbox_inches='tight')

stat.plot_boxplot(np.divide(np.asarray(num_loc_cluster), (0.5*np.asarray(cluster_diameters))**2*np.pi),
                  scale='lin', xlabel='$\,$', bptype='violin',
                  ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
plt.savefig(output_dir + file_name + '_clusterdensities_boxplot.pdf')

