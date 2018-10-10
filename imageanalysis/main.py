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

file_dirs = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/test/pointpattern/synthetic_chiara_Baumgart2016/1/Ch1/']  #,
             #  '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/melike/Histones/crops/TSA/']  # >= 1
is_dataset = 0  # 0 if run analysis with one file ('file_name'); 1 if all files in file_dir(s)
file_name = 'synth_clustered_sig80_mu7_ncl10_r50_Nc50_1'  # if is_dataset=0, then one single file (in file_dirs)
fileExt = '.txt'
is_storm = 0  # 0 if two-columns txt file with xc (x-corrected) and yc (y-corrected) from typical STORM output
run = dict(image_processing=1, image_analysis=0)  # run image_processing and or image_analysis

feature_all = []  # list of all features obtained for each training image in file_dirs
print '\n\n _______START FULL ANALYSIS_______ '; inini_time = time.time()
print ' _________________________________'

if run['image_processing']:
    for file_dir in file_dirs:
        for jj, fileid in enumerate(os.listdir(file_dir)):
            if not is_dataset and (jj > 0): break
            elif is_dataset: file_name = fileid.split(fileExt)[0]

            print '\n\n\n _______ file no.%d _______ ' % jj
            # # ============== INPUT PARAMETERS ========
            # # ========================================
            dict_inputfile = dict(file_dir=file_dir, file_name=file_name, fileExt=fileExt, is_storm=is_storm,
                                                      out_channel='all',
                                  ispp=1, compute_ROI=0, crop=1, crop_range=[85, 120, 105, 135],
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
                data.read(dict_inputfile, save_output_dir=output_dir, plot=True)
                points = data.points
            else:
                image = (cv2.imread(file_dir + file_name + fileExt, 0).astype(float)).T  # recall:

            if dict_inputfile.get('ispp'):
                # util.plot_frameno(data.points, data.frame_no)

                print '\ttotal number of localizations = %d' % points.shape[0]
                vproc.plot_points(points)
                # plt.savefig(output_dir + file_name + '_pp.pdf', bbox_inches='tight')

                # # ============== ROI =============================
                # # ================================================
                print '\n _______PRE-PROCESSING_______'

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
                        vproc.plot_points(points_roi)
                        # plt.savefig(output_dir + file_name + '_pp_roi' + '.pdf', bbox_inches='tight', dpi=400)
                else:
                    print 'No ROI computations required.'
                    points_roi = points  # points_roi = iproc.compute_roi(compute_ROI, points)

                print '\tnumber of localizations in the analysis = %d' % points_roi.shape[0]
                if dict_inputfile['compute_ROI']:
                    print ("\n\n***FINISHED : ROI computed and saved as non-storm file (total time =  %.2f seconds)" % (time.time() - inini_time))
                    exit()

                # # ====== GENERATE IMAGE =====================
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

                    start_time = time.time()
                    print 'Computing Voronoi tessellation... '
                    vor = Voronoi(points_roi)  # compute Voronoi tessellation
                    print time.time()-start_time

                    print 'Computing Voronoi descriptors...'
                    start_time = time.time()
                    vproc.compute_parameters(vor, dict_inputfile)   # new attribute in vor object: vor.areas
                    print time.time()-start_time
                    print "\tnuclear area = %.0f" % float(vor.areas_total), '[o. pixel size]2'

                    start_time = time.time()
                    print 'Converting Voronoi tessellation into image (i.e., interpolate)...',
                    image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                                        interpolate_method=dict_image.get('interpolate_method'),
                                                        fill_value=0.0,
                                                        density_transform=dict_image['detect_densitytransform'])  # log 4 detect
                    print time.time()-start_time
                    image_descr = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                                              interpolate_method=dict_image.get('interpolate_method'),
                                                              fill_value=0.0,
                                                              density_transform=dict_image['descr_densitytransform'])
                    print('Done.')
                    print 'Plotting Voronoi zero-rank densities image...'
                    iproc.plot_image(image, cmap='jet', norm='lin', plot_axis='on')
                    # plt.savefig(output_dir + file_name + '_densities_image.pdf', bbox_inches='tight'); print 'Saved.'

                    # print 'Plotting Voronoi zero-rank densities point pattern...'
                    # threshold = float(2*vor.densities_average); vproc.threshold(vor, thr=threshold)
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

            print 'Intensity-dependent computations...'
            dict_sift = dict(scale_pixel_size=scale_pixel_size, resolution=1,
                             original_pixel_size=dict_inputfile.get('original_pixel_size'),

                             # feature detection
                             t=10, feature_name='blob',  # 0.8, scale_ini=80
                             thresholding=0, threshold_percent=0.5, num_features='all', scale_range_is='nm',
                             scale_ini=50,
                             scale_end=600,  # diam. of search, if pixel->analysispixel
                             scale_spacing='odd', nscales=150,
                             scale_resolution=dict_inputfile.get('resolution'),  # 'scale_resolution': 1, # (radius) in
                             # scale_range_is (if [nm] and STORM -> min. 20nm)
                             max_filter_width=7, max_filter_depth=7,

                             # feature description [main orientation(s)]
                             compute_orientation=1, n_bins_ori=36, peak_ratio=0.8,
                             sigma_ori_times=1.5, window_ori_radtimes=1, smooth_cycles=2,

                             # feature description [histograms]
                             compute_sift_descr=1,
                             sigma_descr_times=2, window_descr_radtimes=1, n_hist=16, n_bins_descr=8, threshold_sat=0.2,
                             plot_graphics=0)
            util.saveparameters(output_dir + file_name + '_parameters.txt', dict1=dict_inputfile, dict2=dict_image, dict3=dict_sift)

            reload(iproc); plt.close('all')
            print 'Detecting and describing features (%s)...' % dict_sift.get('feature_name'); start_time = time.time()
            feature = iproc.find_feature(image, dict_sift, image_descr=image_descr)
            print ("\tDone (total time =  %.2f seconds)" % (time.time() - start_time))
            if dict_inputfile.get('ispp'):
                feature['number_localizations'] = vproc.localizations_feature(vor, feature, dict_sift)
                feature['density'] = feature['number_localizations']/(np.pi*(0.5*feature['diameter'])**2)
            print "\tnumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape

            print 'Plotting intensity-dependent Voronoi features...'
            start_time = time.time()
            iproc.plot_feature(image, feature, cmap='jet', norm='linear', plot_axis='on', blob_color='strength')
            plt.savefig(output_dir + file_name + '_features_image.pdf', bbox_inches='tight')
            print("Done (time =  %.2f seconds)" % (time.time() - start_time))

            if dict_inputfile.get('ispp'):
                print 'Plotting intensity-dependent Voronoi features on Voronoi...',
                start_time = time.time()
                vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on')
                plt.savefig(output_dir + file_name + '_features_pp.pdf', bbox_inches='tight')
                print ("Done (time =  %.2f seconds)" % (time.time() - start_time))

            feature['file_name'] = file_name; feature['image'] = image; feature['vor'] = vor
            feature_all.append(feature)

            dict_output = dict(image_area=image.shape[0]*image.shape[1]*analysis_pixel_size**2,
                               ntotal=points.shape[0],
                               number_clusters=feature.get('argmaxgrad')[0].shape[0],
                               cluster_diameters=feature['diameter'].tolist(),
                               number_localizations_cluster=feature['number_localizations'].tolist(),
                               )
            util.saveparameters(output_dir + file_name + '_output.txt', dict1=dict_output)
            print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

    # # ====== IMAGE PROCESSING - SAMPLE STATISTICS
    # # ===========================================
    print '\n _______SAMPLE STATISTICS_______'; ini_time = time.time()
    reload(vproc); reload(iproc); reload(stat); reload(util)
    print 'Computing feature statistics...'

    if len(feature_all) > 1: savefile_suffix_all = output_dir + '_'.join(file_name.split('_')[0:-1])
    else: savefile_suffix_all = output_dir + file_name
    stat.statistic_descrip(feature_all, ispp=dict_inputfile['ispp'], pixel_size=analysis_pixel_size,
                           savefile_suffix=savefile_suffix_all)
    print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))


# # ====== IMAGE ANALYSIS ====================
# # ==========================================
if run['image_analysis']:
    plt.close('all')
    reload(stat); reload(iproc); reload(vproc); reload(util)
    # # ==== CREATE VOCABULARY, unsupervised =====
    if dict_sift.get('compute_orientation'):
        print '\n _______ IMAGE ANALYSIS_______'; ini_time = time.time()
        init = 'k-means++'
        k0 = 2; kn = 2; sserror = []  # build bag of words - unsupervised kmeans
        for k in range(k0, kn+1):
            print 'Creating vocabulary with %d words...' % k
            histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
            histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in
                                    hist])
            print 'size all histograms: ', len(histogram)
            kmeans = stat.create_vocabulary(feature_all, dict_sift, n_cluster=k, init=init)
            sserror.append(kmeans.inertia_)

            # feature_all_sub = stat.filter_features(feature_all, kmeans, ball_radius=4.5*10**-1)

            # plot the feature class with k+1 according to the colorcode at k. Comment if it feels unnecessary.
            if k == k0: centers_permuted0 = kmeans.cluster_centers_
            if k > k0: centers_permuted0, kmeans = util.permute_labels(kmeans, k, centers_permuted0)

            image_no = 0  # len(feature_all)-1
            feature_image = feature_all[image_no]
            image = feature_image['image']
            vor = feature_image['vor']
            if len(feature_all) > 1:
                savefile_suffix_all = output_dir + '_'.join(file_name.split('_')[0:-1]) + '_kmeans%d' % k
            else:
                savefile_suffix_all = output_dir + feature_image['file_name'] + '_kmeans%d' % k
            stat.scatterplot_vocabulary(feature_all, kmeans, n_cluster=k, cmap=util.discrete_cmap(k, "jet"),
                                        savefile_suffix=savefile_suffix_all)
            stat.sift2shape(kmeans, cmap=util.discrete_cmap(k, "jet"), savefile_suffix=savefile_suffix_all)

            stat.words_characteristics(feature_all, kmeans.labels_, cmap=util.discrete_cmap(k, "jet"),
                                       pixel_size=analysis_pixel_size, savefile_suffix=savefile_suffix_all)

            # plots and statistics on one selected image, optional (change kmeans.labels_ -> kmeans_labels_image):
            kmeans_labels_image = util.select_labels_image(feature_all, kmeans.labels_, image_no=image_no)
            savefile_suffix = output_dir + feature_image['file_name'] + '_kmeans%d' % k  # or None
            stat.words_characteristics(feature_image, kmeans_labels_image, cmap=util.discrete_cmap(k, "jet"),
                                       pixel_size=analysis_pixel_size, savefile_suffix=savefile_suffix)
            iproc.plot_feature(image, feature_image, cmap='gray', norm='linear', plot_axis='on',
                               blob_color='class', ori_color=kmeans_labels_image, ori_cmap=util.discrete_cmap(k, "jet"))
            plt.savefig(savefile_suffix + '_features_image.pdf', bbox_inches='tight')
            if dict_inputfile['ispp']:
                vproc.plot_feature(vor, feature_image, dict_sift, show_points=True, cmap='gray', norm='log',
                                   plot_axis='on', blob_color='class', ori_color=kmeans_labels_image,
                                   ori_cmap=util.discrete_cmap(k, "jet"))
                plt.savefig(savefile_suffix + '_features_pp.pdf', bbox_inches='tight')

        fig2, ax2 = plt.subplots()
        ax2.plot(range(k0, kn+1), sserror, 'k.-', markersize=10)
        ax2.set_ylabel('total within sum of squares'); ax2.set_xlabel(r'number of clusters $k$')
        fig2.savefig('_'.join(savefile_suffix_all.split('_')[0:-1]) + '_withinsumofsquares.pdf',
                     bbox_inches='tight')
        print ("Done (total time =  %.2f seconds)" % (time.time() - ini_time))

print ("\n\n***FINISHED ANALYSIS (total time =  %.2f seconds)" % (time.time() - inini_time))

# # if __name__ == '__main__':
# #     main()
