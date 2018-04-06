# modules
import sys, os, time
import numpy as np, matplotlib.pyplot as plt
import cv2
from scipy.spatial import Voronoi
# my_modules
from src import vprocessing as vproc, iprocessing as iproc, statistics as stat, utilities as util
# from src import *

# def main():
#     # all the stuff

plt.close('all')
# set directories
parent_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis'#; os.chdir(parent_dir)
source_dir = parent_dir + '/src/'; sys.path.insert(0, source_dir)
output_dir = parent_dir + '/output/'  # set output directory
if not os.path.exists(output_dir): os.makedirs(output_dir)

# # =========== READ INPUT FILE =================
# # ================================================
experiment_author = ''; fileDir = ''; fileName = ''

fileDir = '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-18_HeLa_DualColor_SMC1_SMC3/SMC1_SMC3 in DMSO Controls/Ch2/'
fileName = 'SMC3_SMC1_DMSO_000_list_drift_1740-2622_Ch2_145165_105125'

# for file in os.listdir(fileDir):
#     fileName = file.split(".txt")[0]

data = util.pointpattern()
data.read(fileDir, fileName=fileName, fileExt='.txt', storm=0, channels_num=1, out_channel=[1, 2],
                          save_output_dir=None, plot=True)
points = data.points

# image = (cv2.imread('../data/test/myshape00.png', 0).astype(float)).T  # recall: inputfile_ispp = False

# stat.plot_frameno(data.points, data.frame_no)

# # ============== INPUT PARAMETERS ========
# # ========================================
# cell_no = str(0) + str(2)
dict_inputfile = {'filename': fileName,
                  'ispp': 1, 'compute_ROI': 0, 'crop': 1, 'crop_range': [170, 200, 80, 100],
                             'pixelate': 0,
                             'tessellate': 1,
                  'original_pixel_size': 160, 'photonconv': 0.14, 'resolution': 0.1}   # [nm]/[pixel], e.g. STORM
analysis_pixel_size = 20  # [nm] <<< 160 [nm] (STORM res.) -> scale pixel size anal.p.s/ori.p.s
scale_pixel_size = float(analysis_pixel_size)/dict_inputfile.get('original_pixel_size')
dict_image = {'scale_pixel_size': scale_pixel_size,
              'original_pixel_size': dict_inputfile.get('original_pixel_size'),
              'interpolate_method': 'linear'}
# util.saveparameters(output_dir + 'parameters.txt', dict1=dict_inputfile, dict2=dict_image)

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
            roi_scale_pixel_size = 200. / 160  # [pixel] w.r.t original pixel size = [pixel]
            dict_roi = {'scale_pixel_size': roi_scale_pixel_size, 't': 40, 'feature_name': 'edge',
                        'thresholding': False, 'threshold_percent': 0.6}

            image, image_ptslabel = iproc.pattern2image(points, dict_roi.get('scale_pixel_size'))
            image_blurred = iproc.blur_image(image, dict_roi.get('t'))
            feature = iproc.find_feature(image, dict_roi)
            points_roi = iproc.image2pattern(feature.get('image_roi'), points, image_ptslabel)

            # # visualization
            iproc.plot_image(image, cmap='gray', norm='lin'), plt.title(r"pixelated image")
            iproc.plot_image(image_blurred, cmap='gray', norm='lin'), plt.title(r"blurred pixelated image")
            iproc.plot_feature(image, feature, cmap='gray', norm='lin', feature_name=dict_roi.get('feature_name'))
            iproc.plot_image(feature.get('image_roi'), cmap='gray', norm='lin'), plt.title('ROI')
            print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
        else:
            crop_range = dict_inputfile.get('crop_range')
            points_roi = iproc.points_2dcrop(data.points, crop_range)
            # fileName_crop = fileName + '_' + 'Ch1' + '_' + \
            #                 str(crop_range[0])+str(crop_range[1]) + '_' +\
            #                 str(crop_range[2]) + str(crop_range[3])
            # np.savetxt(output_dir + fileName_crop + '.txt', points_roi)
            # points_roi = iproc.points_2dcrop(data.points2, crop_range)
            # fileName_crop = fileName + '_' + 'Ch2' + '_' + \
            #                 str(crop_range[0])+str(crop_range[1]) + '_' +\
            #                 str(crop_range[2]) + str(crop_range[3])
            # np.savetxt(output_dir + fileName_crop + '.txt', points_roi)
    else:
        print '\nNo ROI computations required.'
        points_roi = points  # points_roi = iproc.compute_roi(compute_ROI, points)

    vproc.plot_points(points_roi)

    plt.savefig(output_dir + 'pp_roi' + '.pdf', bbox_inches='tight')
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
        threshold = float(2*vor.densities_average)
        vproc.threshold(vor, thr=threshold)
        vproc.plot_densities(vor, thr=None, show_points=True, cmap='jet', norm='log', plot_axis='on')
        plt.savefig(output_dir + 'densities_pp.pdf', bbox_inches='tight'); print 'Saved.'

        fileName_crop = fileName + '_' + 'nonoise'
        np.savetxt(output_dir + fileName_crop + '.txt', vor.points_thresholded)

        # print 'Plotting Voronoi areas...'
        # threshold = float((2*vor.densities_average)**-1)
        # vproc.plot_areas(vor, threshold=threshold, show_points=False, plot_axis='off')
        # plt.savefig(output_dir + 'areas_pp.pdf', bbox_inches='tight')

        # print '\nPlotting Voronoi tessellation... '
        # dict_plotvoronoi2d = {'show_vertices': 'True', 'show_points': 'True', 'line_width': '0.5',
        #                       'show_unbounded_cells': 'True'}
        # vproc.voronoi_plot_2d(vor, **dict_plotvoronoi2d), plt.title(r"Voronoi tessellation")
        # plt.savefig(output_dir + 'tessellation_pp.pdf', bbox_inches='tight')

        print (" _______VORONOI ANALYSIS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))


# # ====== IMAGE PROCESSING ===================
# # ===========================================
print '\n _______IMAGE PROCESSING_______'; ini_time = time.time()
plt.close('all')
reload(vproc); reload(iproc); reload(stat)

print '\nIntensity-dependent computations...'; start_time = time.time()
dict_sift = {'scale_pixel_size': scale_pixel_size, 'resolution': 1,  # resolution=10??
             'original_pixel_size': dict_inputfile.get('original_pixel_size'),

              # feature detection
             't': 10, 'feature_name': 'blob',
             'thresholding': 1, 'threshold_percent': 0.95,
             'scale_range_is': 'nm', 'scale_ini': 200, 'scale_end': 300, 'scale_spacing': 'odd',  # diameter of search
             'nscales': 150,  # 'scale_resolution': 1, # (radius) in scale_range_is (if [nm] and STORM -> min. 20nm)
                             'scale_resolution': 2*dict_inputfile.get('resolution'),

             'max_filter_width': 3, 'max_filter_depth': 3,  # [scale pixel size]

             # feature description [main orientation(s)]
             'compute_orientation': False, 'n_bins_ori': 36, 'peak_ratio': 0.7, 'sigma_ori_times': 1.5,
             'window_ori_radtimes': 2, 'smooth_cycles': 2,

             # feature description [histograms]
             'compute_sift_descr': True, 'sigma_descr_times': 1.5, 'window_descr_radtimes': 2, 'n_hist': 16,
             'n_bins_descr': 8, 'threshold_sat': 0.2, 'plot_graphics': 0}

util.saveparameters(output_dir + 'parameters.txt', dict1=dict_inputfile, dict2=dict_image, dict3=dict_sift)

print 'Finding intensity-dependent features (%s)...' % dict_sift.get('feature_name')
feature = iproc.find_feature(image, dict_sift)
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
print "\tnumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape

print 'Plotting intensity-dependent Voronoi features... ...'; start_time = time.time()
iproc.plot_feature(image, feature, feature_name=dict_sift.get('feature_name'), cmap='jet', norm='log', plot_axis='on')
plt.savefig(output_dir + 'features_image.pdf', bbox_inches='tight')
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))

print 'Plotting intensity-dependent Voronoi features on Voronoi...'; start_time = time.time()
vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on')
# plt.savefig(output_dir + 'features_pp.pdf', bbox_inches='tight')
print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))

print (" _______IMAGE PROCESSING: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))

# # ====== IMAGE PROCESSING - SAMPLE STATISTICS
# # ===========================================
print '\n _______SAMPLE STATISTICS_______'; ini_time = time.time()
reload(vproc); reload(iproc); reload(stat); reload(util)
print 'Computing intensity-dependent Voronoi features...'
blob_diameters = analysis_pixel_size*3*np.sqrt(feature.get('tnew'))  # in original units
nnd_localizations = iproc.nnd_feature(feature, dict_sift)
if dict_inputfile.get('tessellate'):
    number_localizations = vproc.localizations_feature(vor, feature, dict_sift)
    densities = number_localizations/((0.5*blob_diameters)**2*np.pi)
    vareas_statistics = stat.sample_statistics(vor.areas[vor.areas < float('inf')]*dict_inputfile.get('original_pixel_size')**2)
    feature_statistics = stat.sample_statistics(number_localizations)
    densities_statistics = stat.sample_statistics(densities)
    print '\tVoronoi polygon areas: \t', vareas_statistics, '[nm2]'
    print '\tNumber of loc. per blob: \t', feature_statistics, '[loc/blob]'
    print '\tLocal blob densities: \t', densities_statistics, ' [loc/nm2]'
    stat.plot_hist(number_localizations, hist_scale='log', num_bins=50, xlabel=r'number of localizations per blob')
    stat.plot_hist(vor.areas, hist_scale='log', bins=analysis_pixel_size * 3 * np.sqrt(feature.get('scale_range')),
                   xlabel=r'Voronoi polygon area [nm$^2$]')
    stat.plot_hist(densities, hist_scale='lin', num_bins=50, xlabel=r'blob density [localizations/nm$^2$]')
stat.plot_hist(blob_diameters, bins=analysis_pixel_size*3*np.sqrt(feature.get('scale_range')), xlabel=r'blob diameter [nm]')
stat.errorbar_featureresponse(feature, dict_sift, xlabel=r'blob diameter [nm]')

print (" _______SAMPLE STATISTICS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))

# # ====== IMAGE ANALYSIS ====================
# # ==========================================
print '\n _______ DESCRIPTORS ANALYSIS_______'; ini_time = time.time()
plt.close('all'); reload(stat); reload(iproc)
fig2, ax2 = plt.subplots()
for ii in range(2, 3):
    fig1, ax1 = plt.subplots()
    kmeans = stat.siftdescr_analysis(feature, dict_sift, n_cluster=ii,  max_iter=300, n_jobs=-1,
                                     compute_pca=True, plot_graphics=True, fig=fig1, ax=ax1)
    fig1.savefig(output_dir + 'kmeans%d' % ii + '.pdf', bbox_inches='tight')
    ax2.plot(ii, kmeans.inertia_, 'k.', markersize=10)
    ax2.set_ylabel('total within sum of squares'); ax2.set_xlabel('number of clusters k'); fig2.hold(True)

fig2.savefig(output_dir + 'within_sum_of_squares.pdf', bbox_inches='tight')
iproc.plot_feature(image, feature, feature_name=dict_sift.get('feature_name'), cmap='gray', norm='log', plot_axis='on',
                   cluster_color=kmeans.labels_)
plt.savefig(output_dir + 'features_image_kmeans.pdf', bbox_inches='tight')
plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on',
             cluster_color=kmeans.labels_)

print (" _______DESCRIPTORS ANALYSIS: DONE (total time =  %.2f seconds)" % (time.time() - ini_time))


# if __name__ == '__main__':
#     main()
