# modules
import sys, os, time
import numpy as np, matplotlib.pyplot as plt
import cv2
from scipy.spatial import Voronoi
# my_modules
from src import vprocessing as vproc, iprocessing as iproc, statistics as stat, utilities as util

# def main():
#     # all the stuff

plt.close('all')
# set directories
parent_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis'   # os.chdir(parent_dir)
source_dir = parent_dir + '/src/'; sys.path.insert(0, source_dir)
output_dir = parent_dir + '/output/'  # set output directory
if not os.path.exists(output_dir): os.makedirs(output_dir)

# # =========== READ INPUT FILE =================
# # ================================================
experiment_author = ''; file_dir = ''; file_name = ''

file_dir = '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/test/pointpattern/classification' \
           '/circles_rectangles/synthetic_circlerectangle_enrich20_densityc325px2/'
# file_name = 'synthetic_circlerectangle_enrich20_densityc325px2_2'
fileExt = '.txt'

feature_all = []  # list of all features obtained for each training image in file_dir

for jj, fileid in enumerate(os.listdir(file_dir)):
    # if jj > 0: break
    file_name = fileid.split(fileExt)[0]

    print '\n\n\n _______START ANALYSIS_______';  inini_time = time.time()

    # # ============== INPUT PARAMETERS ========
    # # ========================================
    # cell_no = str(0) + str(2)
    dict_inputfile = dict(file_name=file_name,
                          ispp=1, compute_ROI=0, crop=1, crop_range=[100, 125, 21, 77],
                                  pixelate=0,
                                  tessellate=1,
                          original_pixel_size=160, photonconv=0.14, resolution=1)   # [nm]/[pixel], e.g. STORM
    analysis_pixel_size = 10  # 10 [nm] <<< 160 [nm] (STORM res.) -> scale pixel size anal.p.s/ori.p.s
    scale_pixel_size = float(analysis_pixel_size)/dict_inputfile.get('original_pixel_size')
    dict_image = dict(scale_pixel_size=scale_pixel_size,
                      original_pixel_size=dict_inputfile.get('original_pixel_size'),
                      interpolate_method='linear', detect_densitytransform='log', descr_densitytransform='linear')
    # read
    if dict_inputfile.get('ispp'):
        data = util.pointpattern()
        data.read(file_dir, file_name=file_name, fileExt='.txt', storm=0, channels_num=2, out_channel=[0, 1], plot=True)
        points = data.points
    else:
        image = (cv2.imread(file_dir + file_name + fileExt, 0).astype(float)).T  # recall:

    if dict_inputfile.get('ispp'):

        # util.plot_frameno(data.points, data.frame_no)

        print '\ttotal number of localizations = %d' % points.shape[0]
        # vproc.plot_points(points)
        # plt.savefig(output_dir + file_name + '_pp.pdf', bbox_inches='tight')

        # # ============== ROI =============================
        # # ================================================
        print '\n _______PRE-PROCESSING_______'
        reload(iproc); reload(stat); reload(vproc)

        if dict_inputfile.get('compute_ROI'):
            if not dict_inputfile.get('crop'):
                print 'Computing ROI (scale-space)...'; start_time = time.time()
                roi_scale_pixel_size = 200. / 160  # [pixel] w.r.t original pixel size = [pixel]
                dict_roi = {'scale_pixel_size': roi_scale_pixel_size, 't': 40, 'feature_name': 'edge',
                            'thresholding': False, 'threshold_percent': 0.6}

                image, image_ptslabel = iproc.pattern2image(points, dict_roi.get('scale_pixel_size'))
                image_blurred = iproc.blur_image(image, dict_roi.get('t'))
                feature = iproc.find_feature(image, dict_roi)
                points_roi = iproc.image2pattern(feature.get('image_roi'), points, image_ptslabel)

                # # visualization
                iproc.plot_image(image, cmap='gray', norm='lin'); plt.title(r"pixelated image")
                iproc.plot_image(image_blurred, cmap='gray', norm='lin'), plt.title(r"blurred pixelated image")
                iproc.plot_feature(image, feature, cmap='gray', norm='lin')
                iproc.plot_image(feature.get('image_roi'), cmap='gray', norm='lin'), plt.title('ROI')
                print("\tDONE (time =  %.2f seconds)" % (time.time() - start_time))
            else:
                crop_range = dict_inputfile.get('crop_range')
                points_roi = iproc.points_2dcrop(data.points, crop_range)
                # file_name_crop = file_name + '_'
                # np.savetxt(output_dir + file_name_crop + '.txt', points_roi)
                # points_roi = iproc.points_2dcrop(data.points2, crop_range)
                file_name_crop = file_name + '_' + \
                                str(crop_range[0])+str(crop_range[1]) + '_' +\
                                str(crop_range[2]) + str(crop_range[3])
                np.savetxt(output_dir + file_name_crop + '.txt', points_roi)
        else:
            print 'No ROI computations required.'
            points_roi = points  # points_roi = iproc.compute_roi(compute_ROI, points)

        vproc.plot_points(points_roi)
        # plt.savefig(output_dir + file_name + '_pp_roi' + '.pdf', bbox_inches='tight')
        print '\tnumber of localizations in the analysis = %d' % points_roi.shape[0]

        # # ====== IMAGE GENERATION ===================
        # # ===========================================
        if dict_inputfile.get('pixelate'):
            print 'Pixelating...',
            image = iproc.pattern2image(160*points, pixel_size=analysis_pixel_size)[0]
            iproc.plot_image(image, cmap='jet', norm='log', plot_axis='on')
            # plt.savefig(output_dir + file_name + '_pixelated_pp.pdf', bbox_inches='tight')
            print 'Done.'

        if dict_inputfile.get('tessellate'):
            plt.close('all')
            print '\n _______VORONOI ANALYSIS_______'; ini_time = time.time()
            reload(vproc); reload(iproc)

            print 'Computing Voronoi tessellation... '
            vor = Voronoi(points_roi)  # compute Voronoi tessellation

            print 'Computing Voronoi descriptors...'
            vproc.compute_parameters(vor, dict_inputfile)   # new attribute in vor object: vor.areas
            print "\tnuclear area = %.0f" % float(vor.areas_total), '[o. pixel size]2'

            print 'Converting Voronoi tessellation into image (i.e., interpolate)...',
            image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                                interpolate_method=dict_image.get('interpolate_method'), fill_value=0.0,
                                                density_transform=dict_image['detect_densitytransform'])  # log 4 detect
            image_descr = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                                      interpolate_method=dict_image.get('interpolate_method'),
                                                      fill_value=0.0,
                                                      density_transform=dict_image['descr_densitytransform'])
            print('Done.')
            print 'Plotting Voronoi zero-rank densities image...'
            iproc.plot_image(image, cmap='jet', norm='lin', plot_axis='on')
            # plt.savefig(output_dir + file_name + '_densities_image.pdf', bbox_inches='tight'); print 'Saved.'

            print 'Plotting Voronoi zero-rank densities point pattern...'
            threshold = float(2*vor.densities_average); vproc.threshold(vor, thr=threshold)
            # vproc.plot_densities(vor, thr=None, show_points=True, cmap='jet', norm='log', plot_axis='on')
            # plt.savefig(output_dir + file_name + '_densities_pp.pdf', bbox_inches='tight'); print 'Saved.'

            # print 'Plotting Voronoi areas...'
            # threshold = float((2*vor.densities_average)**-1)
            # vproc.plot_areas(vor, threshold=threshold, show_points=False, plot_axis='off')
            # plt.savefig(output_dir + file_name + '_areas_pp.pdf', bbox_inches='tight')

            # print '\nPlotting Voronoi tessellation... '
            # dict_plotvoronoi2d = {'show_vertices': 'True', 'show_points': 'True', 'line_width': '0.5',
            #                       'show_unbounded_cells': 'True'}
            # vproc.voronoi_plot_2d(vor, **dict_plotvoronoi2d), plt.title(r"Voronoi tessellation")
            # plt.savefig(output_dir + file_name + '_tessellation_pp.pdf', bbox_inches='tight')

            print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

    # # ====== IMAGE PROCESSING ===================
    # # ===========================================
    print '\n _______IMAGE PROCESSING_______'; ini_time = time.time()
    plt.close('all')
    reload(vproc); reload(iproc); reload(stat)

    print 'Intensity-dependent computations...'; start_time = time.time()
    dict_sift = dict(scale_pixel_size=scale_pixel_size, resolution=1,
                     original_pixel_size=dict_inputfile.get('original_pixel_size'),

                     # feature detection
                     t=10, feature_name='blob',  # 0.8, scale_ini=80
                     thresholding=1, threshold_percent=0.6, num_features=250, scale_range_is='nm', scale_ini=70,
                     scale_end=200,  # diam. of search, if pixel->analysispixel
                     scale_spacing='odd', nscales=150,
                     scale_resolution=dict_inputfile.get('resolution'),  # 'scale_resolution': 1, # (radius) in
                     # scale_range_is (if [nm] and STORM -> min. 20nm)
                     max_filter_width=7, max_filter_depth=7,

                     # feature description [main orientation(s)]
                     compute_orientation=1, n_bins_ori=36, peak_ratio=0.7,
                     sigma_ori_times=1.5, window_ori_radtimes=1, smooth_cycles=2,

                     # feature description [histograms]
                     compute_sift_descr=1,
                     sigma_descr_times=2, window_descr_radtimes=1, n_hist=16, n_bins_descr=8, threshold_sat=0.2,
                     plot_graphics=0)
    util.saveparameters(output_dir + file_name + '_parameters.txt', dict1=dict_inputfile, dict2=dict_image, dict3=dict_sift)

    reload(iproc); plt.close('all')
    print 'Finding intensity-dependent features (%s)...' % dict_sift.get('feature_name')
    feature = iproc.find_feature(image, dict_sift, image_descr=image_descr)
    print "\tnumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape
    print("Done (time =  %.2f seconds)" % (time.time() - start_time))

    print 'Plotting intensity-dependent Voronoi features...'
    start_time = time.time()
    iproc.plot_feature(image, feature, cmap='jet', norm='lin', plot_axis='on', blob_color='strength')
    plt.savefig(output_dir + file_name + '_features_image.pdf', bbox_inches='tight')
    print("Done (time =  %.2f seconds)" % (time.time() - start_time))

    if dict_inputfile.get('ispp'):
        print 'Plotting intensity-dependent Voronoi features on Voronoi...',
        start_time = time.time()
        vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on')
        # plt.savefig(output_dir + file_name + '_features_pp.pdf', bbox_inches='tight')
        print ("Done (time =  %.2f seconds)" % (time.time() - start_time))

    print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

    # # ====== IMAGE PROCESSING - SAMPLE STATISTICS
    # # ===========================================
    print '\n _______SAMPLE STATISTICS_______'; ini_time = time.time()
    reload(vproc); reload(iproc); reload(stat); reload(util)
    print 'Computing feature statistics...'
    blob_diameters = analysis_pixel_size*3*np.sqrt(feature.get('scale'))  # in original units
    # nnd_localizations = iproc.nnd_feature(feature, dict_sift)
    stat.plot_hist(0.5*blob_diameters, num_bins=np.unique(feature.get('scale')).shape[0]+1, xlabel=r'blob radius R [nm]')
    # plt.savefig(output_dir + file_name + '_blobradius_hist.pdf', bbox_inches='tight')
    stat.plot_boxplot(0.5*blob_diameters, bptype='violin', ylabel=r'blob radius R [nm]')
    # plt.savefig(output_dir + file_name + '_blobradius_boxplot.pdf', bbox_inches='tight')
    # stat.plot_hist(np.pi*(blob_diameters/2.)**2, num_bins=50, xlabel=r'blob area [nm2]')
    stat.plot_boxplot(np.pi*(blob_diameters/2.)**2, bptype='violin', ylabel=r'blob area [nm2]')
    # plt.savefig(output_dir + file_name + '_blobareas_boxplot.pdf', bbox_inches='tight')
    # stat.errorbar_featureresponse(feature, dict_sift, xlabel=r'blob diameter [nm]')
    print '\tNumber of clusters:\t', feature.get('argmaxgrad')[0].shape[0]
    print '\tDensity of clusters:\t', feature.get('argmaxgrad')[0].shape[0]/float(analysis_pixel_size*image.shape[
        0]*image.shape[1]), '[cluster/nm2]'
    if dict_inputfile.get('ispp'):
        feature['number_localizations'] = vproc.localizations_feature(vor, feature, dict_sift)
        densities = feature['number_localizations']/((0.5*blob_diameters)**2*np.pi)
        percentage_number_localizations = feature['number_localizations']/points.shape[0]*100
        vareas_statistics = stat.sample_statistics(vor.areas[vor.areas < float('inf')]*dict_inputfile.get('original_pixel_size')**2)
        densities_statistics = stat.sample_statistics(densities)
        feature_statistics = stat.sample_statistics(feature['number_localizations'])
        print '\tPercentage of localizations in clusters:\t', np.sum(percentage_number_localizations), '%'
        print '\tVoronoi polygon areas: \t', vareas_statistics, '[nm2]'
        print '\tNumber of loc. per blob: \t', feature_statistics, '[loc/blob]'
        print '\tLocal blob densities: \t', densities_statistics, ' [loc/nm2]'
        # stat.plot_hist(vor.areas, hist_scale='log', num_bins=50, xlabel=r'Voronoi polygon area [nm$^2$]')
        stat.plot_boxplot(vor.areas, scale='log', bptype='violin', ylabel=r'Voronoi polygon area [nm$^2$]')
        # plt.savefig(output_dir + file_name + '_voronoiareas_boxplot.pdf', bbox_inches='tight')
        # stat.plot_hist(number_localizations, hist_scale='lin', num_bins=50, xlabel=r'number of localizations per blob')
        stat.plot_boxplot(number_localizations, scale='lin', bptype='violin', ylabel=r'N$^{cluster}$ [points/cluster]')
        # plt.savefig(output_dir + file_name + '_blobnumloc_boxplot.pdf', bbox_inches='tight')
        # stat.plot_hist(densities, hist_scale='lin', num_bins=50, xlabel=r'blob density [localizations/nm$^2$]')
        stat.plot_boxplot(densities, bptype='violin', ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
        # plt.savefig(output_dir + file_name + '_blobdensities_boxplot.pdf', bbox_inches='tight')

    dict_output = dict(image_area=image.shape[0]*image.shape[1]*analysis_pixel_size**2,
                       ntotal=points.shape[0],
                       number_clusters=feature.get('argmaxgrad')[0].shape[0],
                       cluster_diameters=blob_diameters.tolist(),
                       number_localizations_cluster=number_localizations.tolist(),
                       )
    util.saveparameters(output_dir + file_name + '_output.txt', dict1=dict_output)
    print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

    feature_all.append(feature)
    # plt.close('all')

plt.close('all')
# # ====== IMAGE ANALYSIS ====================
# # ==========================================

reload(stat); reload(iproc); reload(vproc); reload(util)
# # ==== CREATE VOCABULARY, unsupervised =====
if dict_sift.get('compute_orientation'):
    print '\n _______ IMAGE ANALYSIS_______'; ini_time = time.time()
    init = 'k-means++'
    k0 = 2; kn = 3; sserror = []  # build bag of words - unsupervised kmeans
    for k in range(k0, kn+1):
        print '\nCreating vocabulary with %d words...' % k
        kmeans = stat.create_vocabulary(feature_all, dict_sift, n_cluster=k, init=init)
        sserror.append(kmeans.inertia_)

        # plot the feature class with k+1 according to the colorcode at k. Comment if it feels unnecessary.
        if k == k0: centers_permuted0 = kmeans.cluster_centers_
        if k > k0: centers_permuted0, kmeans = util.permute_labels(kmeans, k, centers_permuted0)

        fig1, ax1 = plt.subplots()
        stat.scatterplot_vocabulary(feature_all, kmeans, n_cluster=k, fig=fig1, ax=ax1,
                                    cluster_cmap=util.discrete_cmap(k, "jet"))
        fig1.savefig(output_dir + file_name + '_kmeans%d' % k + '.pdf', bbox_inches='tight')

        # plots and statistics on one selected image, optional (change kmeans.labels_ -> kmeans_labels_image):
        kmeans_labels_image = util.select_labels_image(feature_all, kmeans.labels_, image_no=len(feature_all)-1)
        stat.words_characteristics(feature_all[-1], kmeans_labels_image, cmap=util.discrete_cmap(k, "jet"),
                                   pixel_size=analysis_pixel_size,
                                   savefile_suffix=output_dir + file_name + '_kmeans%d' % k)
        iproc.plot_feature(image, feature_all[-1], cmap='gray', norm='linear', plot_axis='on',
                           blob_color='class', ori_color=kmeans_labels_image, ori_cmap=util.discrete_cmap(k, "jet"))
        plt.savefig(output_dir + file_name + '_kmeans%d' % k + '_features_image.pdf', bbox_inches='tight')
        if dict_inputfile['ispp']:
            vproc.plot_feature(vor, feature_all[-1], dict_sift, show_points=True, cmap='gray', norm='log',
                               plot_axis='on', blob_color='class', ori_color=kmeans_labels_image,
                               ori_cmap=util.discrete_cmap(k, "jet"))
            plt.savefig(output_dir + file_name + '_kmeans%d' % k + '_features_pp.pdf', bbox_inches='tight')

    fig2, ax2 = plt.subplots()
    ax2.plot(range(k0, kn+1), sserror, 'k.-', markersize=10)
    ax2.set_ylabel('total within sum of squares'); ax2.set_xlabel(r'number of clusters $k$')
    fig2.savefig(output_dir + file_name + '_within_sum_of_squares.pdf', bbox_inches='tight')
    print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

print ("\n\n***FINISHED ANALYSIS (total time =  %.2f seconds)" % (time.time() - inini_time))

# # if __name__ == '__main__':
# #     main()
