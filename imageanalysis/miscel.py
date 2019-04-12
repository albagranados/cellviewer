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

# means boxplot reduce randomly points
numloc = [0.7542871445147704, 0.7346207257906975, 0.7141239121992324, 0.6364731033894677, 0.6214606268402468, 0.5543857876123826, 0.5088876915125657, 0.4529232705786746]
area = [2.7360167752522626, 2.7764709965322503, 2.7751240048295736, 2.723448409931618, 2.7298550204200445,
  2.771271994965832, 2.754363434116574, 2.8588060144199914]
strength = [0.37510327019165723, 0.40886251146370167, 0.39428130566030944, 0.36952028614279336, 0.36806713362986504, 0.325915702008907, 0.3181416693462508, 0.2877525010441507]
nnd = [1.0469297237932893, 1.0054832554838478, 1.0330510885494062, 1.040892540685331, 1.0547337860161599, 1.0659529867467021, 1.0834818495718617, 1.1238019657469287]
blobsinblobs = [0.31241934143226485, 0.3063072200025862, 0.3569023016625682, 0.32555140952275397,
                0.31378558144085256, 0.3129958604247032, 0.29982017377606807, 0.285980050335271]
blobdensity = [-2.1146194535722174, -2.1271073048158216, -2.178703239744569, -2.2081024603160526, -2.2542452671552575, -2.3977701217055163, -2.4311842691493175, -2.686701891341101]
blobdensitykappa = [-3.6964327362706575, -3.7815574990755545, -3.7714178677246752, -3.771417867724675, -3.81042294056329, -3.8240344230409566, -3.9076376147426637, -4.059286417330685]


import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams["text.usetex"] = True; mpl.rcParams['font.family'] = 'serif'  # configure latex plots
mpl.rcParams.update({'font.size': 24})
fig, ax = plt.subplots()
plt.plot(area, color='k', linestyle='-',  marker='^', label=r'blob area'); plt.hold(1)
plt.plot(blobsinblobs, '--ko', label=r'blobs in blobs')
plt.plot(blobdensity, color='k', linestyle='-.', marker='v',  label=r'blob density')
plt.plot(nnd, color='k', linestyle=':',  marker='s', label=r'nnd')
plt.plot(numloc, color='k',  marker='<', label=r'points per blob')
plt.plot(blobdensitykappa, color='k', marker='>',  label=r'blobs per area')
ax.legend()
xticklabels = [r'10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%']
ylabel = r'[$10**$]'
plt.setp(ax, xticks=[y for y in range(len(area))], xticklabels=xticklabels, ylabel=ylabel)
plt.savefig(savefile_suffix + '_reduced.pdf', bbox_inches='tight')





