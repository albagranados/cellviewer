# Date: 12 / 2016

# modules
import sys, os, time
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
import cv2  # RUN WITH ASTYPE(FLOAT)
# configure latex plots
matplotlib.rcParams['text.usetex'] = True; matplotlib.rcParams['font.family'] = 'serif'

# module info
__author__ = 'alba granados'; __email__ = "alba.granados@crg.eu"; __project__ = 'voronoi analysis'

# set working directory
parent_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis'; os.chdir(parent_dir)
source_dir = parent_dir + '/src/'; sys.path.insert(0, source_dir)
output_dir = parent_dir + '/output/'  # set output directory
if not os.path.exists(output_dir): os.makedirs(output_dir)

# my_modules
import vprocessing as vproc
import iprocessing as iproc  # image processings
import statistics as stat
import utilities as util

# # =========== READ INPUT FILE =================
# # ================================================
experiment_author = ''; fileDir = ''; fileName = ''

# experiment_author = 'petra'
# fileDir = '/home/alba/Dropbox (CRG)/postdoc_CRG/coding/cellviewer/data/petra/'
# fileName = '110716_c1_NB_pl3_84.8_data_to20000_small'

# fileDir = '/home/alba/Dropbox (CRG)/postdoc_CRG/coding/cellviewer/data/test/synthetic_pp/'
# fileName = 'pp_triangle'


fileDir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/pabloguillaume/'
fileName = 'Cell2_000_no_bg_indiv_multi_SVC_list_Drift_corrected_crop_AND_lyso_list_dc_woch9_crop2'
points = util.read_points(fileDir, fileName=fileName, fileExt='.txt', storm=0, channels_num=2, out_channel='all',
                          save_output_dir=None, plot=True)


# experiment_author = 'victoria'
# fileDir = '/home/alba/Dropbox (CRG)/Shared files Victoria-Alba/2017-12-14_HeLa_DualColor_RNApolII_SMC1/RNApolII_SMC1 in HeLa ActD treated/'
# points = util.read_points(fileDir, fileName=None, fileExt='.txt', storm=True, channels_num=2, save_output_dir=output_dir,
#                          plot=True)

# image = cv2.imread('../data/test/pixelated_triangle_square', 0).astype(float) # inputfile_ispp = False
# image = image.T


# # ============== INPUT PARAMETERS ========
# # ========================================
dict_inputfile = {'filename': fileName,
                  'ispp': 1, 'compute_ROI': 0, 'crop': 1, 'crop_range': [0, 30, 85, 125],
                             'pixelate': 1,
                             'tessellate': 0,
                  'original_pixel_size': 1, 'photonconv': 0.14}   # [nm]/[pixel], e.g. STORM
analysis_pixel_size = 0.2  # [nm] <<< 160 [nm] (STORM res.) -> scale pixel size anal.p.s/ori.p.s
scale_pixel_size = float(analysis_pixel_size)/dict_inputfile.get('original_pixel_size')
dict_image = {'scale_pixel_size': scale_pixel_size, 'resolution': 1, # resolution=10??
              'original_pixel_size': dict_inputfile.get('original_pixel_size'),
              'interpolate_method': 'linear'}
util.saveparameters(output_dir + 'parameters.txt', dict1=dict_inputfile, dict2=dict_image)

if dict_inputfile.get('ispp'):

    print '\ttotal number of localizations = %d\n' % points.shape[0]

    # vproc.plot_points(points)
    # plt.savefig(output_dir + 'pp.pdf', bbox_inches='tight')

    # # ============== ROI =============================
    # # ================================================
    print '\n _______PRE-PROCESSING_______'
    reload(iproc); reload(stat); reload(vproc)

    if dict_inputfile.get('compute_ROI'):
        if not dict_inputfile.get('crop'):
            print '\nComputing ROI (scale-space)...'; start_time = time.time()
            roi_scale_pixel_size = 80. / 160  # [pixel] w.r.t original pixel size = [pixel]
            dict_roi = {'scale_pixel_size': roi_scale_pixel_size, 't': 100, 'feature_name': 'edge',
                        'thresholding': False, 'threshold_percent': 0.6}

            image, image_ptslabel = iproc.pattern2image(points, dict_roi.get('scale_pixel_size'))
            image_blurred = iproc.blur_image(image, dict_roi.get('t'))  # filter image & plot
            feature = iproc.find_feature(image, dict_roi)   # detect feature for ROI extraction & plot
            points_roi = iproc.image2pattern(feature.get('image_roi'), points, image_ptslabel)

            # # visualization
            iproc.plot_image(image), plt.title(r"pixelated image")  # plot low resolution image
            # iproc.plot_image(image_blurred, plot_axis='off'), plt.title(r"blurred pixelated image") # plot blurred image
            iproc.plot_feature(image, feature, dict_roi.get('feature_name'))  # plot detected feature (low-reso image)
            # iproc.plot_image(feature.get('image_roi'))  # plot ROI (image)
            print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
        else:
            points_roi = iproc.points_2dcrop(points, dict_inputfile.get('crop_range'))
            np.savetxt(output_dir + fileName + '_crop.txt', points_roi)
    else:
        print '\nNo ROI computations required.'
        points_roi = points  # points_roi = iproc.compute_roi(compute_ROI, points)

    vproc.plot_points(points_roi)
    plt.savefig(output_dir + 'pp_roi.pdf', bbox_inches='tight')
    print '\tnumber of localizations in the analysis = %d\n' % points_roi.shape[0]

    # # ====== IMAGE GENERATION ===================
    # # ===========================================
    if dict_inputfile.get('pixelate'):
        print 'Pixelating...',
        image = iproc.pattern2image(points, pixel_size=analysis_pixel_size)[0]
        iproc.plot_image(image, cmap='jet', norm='log', plot_axis='on')
        plt.savefig(output_dir + 'pixelated_pp.pdf', bbox_inches='tight')
        print 'Done.'

    if dict_inputfile.get('tessellate'):
        plt.close('all')
        print '\n _______VORONOI ANALYSIS_______'; ini_time = time.time()
        reload(vproc); reload(iproc)

        print 'scale pixel size for intensity-dependent Voronoi analysis = ', scale_pixel_size

        print '\nComputing Voronoi tessellation... '; start_time = time.time()
        vor = Voronoi(points_roi)  # compute Voronoi tessellation
        print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))

        print '\nComputing Voronoi descriptors...'; start_time = time.time()
        vproc.compute_parameters(vor, dict_inputfile)   # new attribute in vor object: vor.areas
        print "\tnuclear area = %.0f" % float(vor.areas_total), '[o. pixel size]2'
        print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))

        print 'Converting Voronoi tessellation into image (i.e., interpolate)...',
        image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                            interpolate_method=dict_image.get('interpolate_method'), fill_value=0.0)
        print 'Plotting Voronoi zero-rank densities image...',
        iproc.plot_image(image, cmap='jet', norm='log', plot_axis='on')
        plt.savefig(output_dir + 'densities_image.pdf', bbox_inches='tight'); print 'Saved.'

        print 'Plotting Voronoi zero-rank densities point pattern...',
        threshold = None  # float(2*vor.densities_average)
        vproc.plot_densities(vor, threshold=threshold, show_points=True, cmap='jet', norm='log', plot_axis='on')
        # plt.savefig(output_dir + 'densities_pp.pdf', bbox_inches='tight'); print 'Saved.'

        # print 'Plotting Voronoi areas...'
        # threshold = float((2*vor.densities_average)**-1)
        # vproc.plot_areas(vor, threshold=threshold, show_points=False, plot_axis='off')
        # plt.savefig(output_dir + 'areas_pp.pdf', bbox_inches='tight')

        # print '\nPlotting Voronoi tessellation... '
        # dict_plotvoronoi2d = {'show_vertices': 'True', 'show_points': 'True', 'line_width': '0.5',
        #                       'show_unbounded_cells': 'True'}
        # vproc.voronoi_plot_2d(vor, **dict_plotvoronoi2d), plt.title(r"Voronoi tessellation")
        # plt.savefig(output_dir + 'tessellation_pp.pdf', bbox_inches='tight')

        # print 'Plotting Voronoi area histogram...'
        # bins = analysis_pixel_size*3*np.sqrt(feature.get('scale_range'))
        # stat.plot_hist(vor.areas, hist_scale='log', num_bins=50); plt.xlabel(r'Voronoi polygon area [nm$^2$]')

        print ("\n _______VORONOI ANALYSIS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))


# # ====== IMAGE PROCESSING ===================
# # ===========================================
print '\n _______IMAGE PROCESSING_______'; ini_time = time.time()
plt.close('all')
reload(vproc); reload(iproc)

print '\nIntensity-dependent computations...'; start_time = time.time()
dict_sift = {'scale_pixel_size': scale_pixel_size, 'resolution': 1, # resolution=10??
             'original_pixel_size': dict_inputfile.get('original_pixel_size'),

              # feature detection
             't': 10, 'feature_name': 'blob',
             'thresholding': True, 'threshold_percent': 0.99,
             'scale_range_is': 'pixel', 'scale_ini': 10, 'scale_end': 40, 'scale_spacing': 'lin',
             'automatic_nscales': 8,

             'max_filter_width': 9, 'max_filter_depth': 7, # [scale pixel size]

             # feature description [main orientation(s)]
             'compute_orientation': True, 'n_bins_ori': 36, 'peak_ratio': 0.7, 'sigma_ori_times': 1.5,
             'window_ori_radtimes': 2, 'smooth_cycles': 2,

             # feature description [histograms]
             'compute_sift_descr': True, 'sigma_descr_times': 1.5, 'window_descr_radtimes': 2, 'n_hist': 16,
             'n_bins_descr': 8, 'threshold_sat': 0.2, 'plot_graphics': False}

util.saveparameters(output_dir + 'parameters.txt', dict1=dict_inputfile, dict2=dict_image, dict3=dict_sift)

print 'Finding intensity-dependent features (%s)...' % dict_sift.get('feature_name')
reload(vproc); reload(iproc); reload(stat)
feature = iproc.find_feature(image, dict_sift)
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
print "\tnumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape

print 'Plotting intensity-dependent Voronoi features... ...'; start_time = time.time()
iproc.plot_feature(image, feature, feature_name=dict_sift.get('feature_name'), cmap='jet', norm='log', plot_axis='on')
plt.savefig(output_dir + 'features_image.pdf', bbox_inches='tight')
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))

reload(vproc)
print 'Plotting intensity-dependent Voronoi features on Voronoi...'; start_time = time.time()
vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on')
# plt.savefig(output_dir + 'features_pp.pdf', bbox_inches='tight')
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))


# # ====== IMAGE PROCESSING - SAMPLE STATISTICS
# # ===========================================

# print 'sample statistics:'
# print '\nCalculating number of localizations per intensity-dependent (blob) cluster...'; start_time = time.time()
# number_localizations = vproc.localizations_feature(vor, feature, dict_sift)
# print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
#
# print '\nCalculating nnd for intensity-dependent cluster...'; start_time = time.time()
# nnd_localizations = vproc.nnd_feature(feature, dict_image)
# print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
#
# blob_diameters = analysis_pixel_size*3*np.sqrt(feature.get('tnew'))
# densities = number_localizations/((0.5*blob_diameters)**2*np.pi)
#
# vareas_statistics = stat.sample_statistics(vor.areas*dict_inputfile.get('original_pixel_size')**2)
# feature_statistics = stat.sample_statistics(number_localizations)
# densities_statistics = stat.sample_statistics(10**3*densities)
# print 'Voronoi polygon areas [nm2]: \t', vareas_statistics
# print 'number of loc. per (blob) feature: \t', feature_statistics
# print 'local (blob) feature densities: \t', densities_statistics, ' [*10(-3) loc/nm2]'
# # print 'Synthetic data info: \n', dict_synthetic_data

# print 'Plotting histograms of intensity-dependent Voronoi features...'
# stat.plot_hist(blob_diameters, bins=analysis_pixel_size*3*np.sqrt(feature.get('scale_range')))
# plt.xlabel(r'blob diameter [nm]')
# stat.plot_hist(feature.get('tnew'), bins=feature.get('scale_range'))
# plt.xlabel(r'scale t [=$\sigma^2$]')

# plt.figure(); plt.hist(blob_diameters, bins=analysis_pixel_size*3*np.sqrt(feature.get(
#     'scale_range')))
# plt.figure(); plt.hist(feature.get('tnew'), bins=feature.get('scale_range'))

# reload(stat)
# print 'Plotting errorbar blob response vs blob diameter...'
# stat.errorbar_featureresponse(feature, dict_analysis)

# print 'Plotting histograms of number of localizations per intensity-dependent Voronoi feature...'
# stat.plot_hist(number_localizations, hist_scale='log', num_bins=50)  # number_localizations.shape[0])
# plt.xlabel(r'number of localizations per (blob) cluster')
# plt.ylim(0, 0.18)

# print 'Plotting histograms of intensity-dependent Voronoi feature densities...'
# stat.plot_hist(densities, hist_scale='lin', num_bins=50)
# plt.xlabel(r'cluster (blob) density [localizations/nm$^2$]')

# plt.figure()
# plt.hist(densities, bins='rice', histtype='step',  color='k'); plt.ylabel(r'counts')
# plt.xlabel(r'cluster (blob) density [localizations/nm$^2$]'); plt.xlim(0, 0.04)

print ("\n _______IMAGE PROCESSING: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))


# # ====== IMAGE ANALYSIS ====================
# # ==========================================
plt.close('all'); reload(stat); reload(iproc)
print '\n _______ DESCRIPTORS ANALYSIS_______'; ini_time = time.time()
fig2, ax2 = plt.subplots()
for ii in range(3, 4):
    fig1, ax1 = plt.subplots()
    kmeans = stat.siftdescr_analysis(feature, dict_sift, n_cluster=ii,  max_iter=300, n_jobs=-1,
                                     compute_pca=True, plot_graphics=True, fig=fig1, ax=ax1)
    fig1.savefig(output_dir + 'kmeans%d'%ii + '.pdf', bbox_inches='tight')
    ax2.plot(ii, kmeans.inertia_, 'k.', markersize=10)
    ax2.set_ylabel('total within sum of squares'); ax2.set_xlabel('number of clusters k'); fig2.hold(True)

fig2.savefig(output_dir + 'within_sum_of_squares.pdf', bbox_inches='tight')
iproc.plot_feature(image, feature, feature_name=dict_sift.get('feature_name'), cmap='gray', norm='log', plot_axis='on',
                   cluster_color=kmeans.labels_)
plt.savefig(output_dir + 'features_image_kmeans.pdf', bbox_inches='tight')
vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on',
                   cluster_color=kmeans.labels_)

print ("\n _______SIFT ANALYSIS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))


# # ========== TESTS with synthetic images =================
# # ========================================================
image = cv2.imread('../data/test/sunflower.jpg', 0).astype(float)

plt.close('all'); reload(util)
plt.figure(); plt.imshow(image.T, interpolation='none', cmap='gray', origin='lower')

dict_inputfile = {'original_pixel_size': 1}   # nm/pixel
dict_analysis = {'analysis_pixel_size': 1}
scale_pixel_size = float(dict_analysis.get('analysis_pixel_size'))/dict_inputfile.get('original_pixel_size')

plt.close('all')
reload(iproc)
dict_sift = {'scale_pixel_size': scale_pixel_size, 'resolution': 1, # resolution=10??
             'original_pixel_size': dict_inputfile.get('original_pixel_size'),
             't': 10, 'feature_name': 'blob',
             'thresholding': False, 'threshold_percent': 0.9, 'interpolate_method': 'cubic',
             'scale_ini': 2, 'scale_end': 10, 'automatic_nscales': 150,
             'scale_spacing': 'odd', 'scale_range_is': 'pixel', 'max_filter_width': 9,
             'max_filter_depth': 100,
             'compute_orientation': False, 'n_bins_ori': 36, 'peak_ratio': 0.9, 'sigma_ori_times': 1.5,
             'window_ori_radtimes': 2, 'smooth_cycles': 2,
             'compute_sift_descr': True, 'sigma_descr_times': 1.5, 'window_descr_radtimes': 2, 'n_hist': 4,
             'n_bins_descr': 8, 'threshold_sat': 0.2,
             'plot_graphics': False}
feature = iproc.find_feature(image, dict_sift)
print "\nNumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape
iproc.plot_feature(image, feature, dict_sift.get('feature_name'), cmap='gray', norm='lin',
                   blob_color='r')

reload(stat)
bins = dict_analysis.get('analysis_pixel_size')*3*np.sqrt(feature.get('scale_range'))
blob_diameters = dict_analysis.get('analysis_pixel_size')*3*np.sqrt(feature.get('tnew'))
stat.plot_hist(blob_diameters, bins=np.append(bins, bins[-1]+3))
plt.xlabel(r'blob diameter [pix]')

# average_strength = stat.errorbar_featureresponse(feature, dict_analysis)
# hist, bin_edges = np.histogram(dict_analysis.get('analysis_pixel_size')*3*np.sqrt(feature.get('tnew')),
#                                bins=np.append(dict_analysis.get('analysis_pixel_size')*3*np.sqrt(feature.get(
#                                'scale_range')),feature.get('scale_range')[-1]+0.001))
# counts_norm = np.multiply(hist, average_strength)/np.sum(average_strength)
# plt.figure(); plt.plot(dict_analysis.get('analysis_pixel_size')*3*np.sqrt(feature.get('scale_range')), counts_norm,
#  'k.', markersize=10)

# stat.plot_hist(feature.get('tnew'), bins=feature.get('scale_range'))
# plt.xlabel(r'scale t [=$\sigma^2$]')
# plt.figure(); plt.hist(feature.get('tnew'), bins=feature.get('scale_range'))

# plt.imsave(output_dir + 'pp_triangle_square.png', image)
# # ---------------------------------------

# # ====== DESCRIPTORS ANALYSIS =================
# # =============================================
plt.close('all')
reload(stat), reload(iproc)
print '\n _______ DESCRIPTORS ANALYSIS_______'; ini_time = time.time()
# plt.figure();
for ii in range(2, 4):
    kmeans = stat.siftdescr_analysis(feature, dict_sift, n_cluster=ii,  max_iter=300, n_jobs=-1, compute_pca=True)
    # plt.plot(ii, kmeans.inertia_, 'k+', markersize=10); plt.title('sum of squared distances to the closest centroid')
    # plt.ylabel('total within sum of squares'); plt.xlabel('number of clusters k')
    # plt.hold(True)
iproc.plot_feature(image, feature, dict_sift.get('feature_name'), cmap='gray', norm='lin',
                   cluster_color=kmeans.labels_)  # blob_color='r',

print ("\n _______SIFT ANALYSIS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))
