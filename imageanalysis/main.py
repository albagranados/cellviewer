# modules
import sys, os, time
import numpy as np, matplotlib.pyplot as plt
import cv2, pickle
from scipy.spatial import Voronoi
import matplotlib as mpl
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

# exp_names = ['waplkd_dmso', 'waplkd_actd']
# file_dirs = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/vicky/2018-05-18 aSMC1Rb HeLa Cas9 WAPLKO DMSO ActD/aSMC1 WAPL KO ActD AND aSMC1 WAPL KO DMSO_WINDOWS/aSMC1 WAPL KO DMSO/',
#              '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/vicky/2018-05-18 aSMC1Rb HeLa Cas9 WAPLKO DMSO ActD/aSMC1 WAPL KO ActD AND aSMC1 WAPL KO DMSO_WINDOWS/aSMC1 WAPL KO ActD/',
#              ]
exp_names = ['mitoc', 'microt', 'lyso']
file_dirs = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/pablo/mithocondria/',
             '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/pablo/microtubules/',
             '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/data/guillaume/moleucles list of '
             'lysosomes/roi/']
is_dataset = 1  # 0 if run analysis with one file ('file_name'); 1 if all files in fil'e_dir(s)
file_name = '1_lyso_006_lis_dct_woch9_wobeads'
fileExt = '.txt'
is_storm = 0  # 0 if two-columns txt file with xc (x-corrected) and yc (y-corrected) from typical STORM output
run = dict(image_processing=0, image_analysis=1)  # run image_processing and or image_analysis

if not run['image_processing']:
    reload(util)
    input_files = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/mitochondria'
                    '/mithoc_variables',
                    '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/microtubules'
                    '/microt_variables',
                    '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/lysosome'
                    '/lyso_variables'
                   ]
    # input_files = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/test/pointpattern/circle_rectangle_lownoise/circles/circles_variables',
    #                '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/test/pointpattern/circle_rectangle_lownoise/rectangles/rectangles_variables']
    # input_files = ['/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/vicky/2018-05-18 '
    #                'aSMC1Rb HeLa Cas9 WAPLKO DMSO ActD/aSMC1 WAPL KO ActD AND aSMC1 WAPL KO DMSO_WINDOWS/aSMC1 WAPL KO DMSO/waplkd_dmso_variables',
    #                '/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/imageanalysis/output/vicky/2018-05-18 aSMC1Rb HeLa Cas9 WAPLKO DMSO ActD/aSMC1 WAPL KO ActD AND aSMC1 WAPL KO DMSO_WINDOWS/aSMC1 WAPL KO ActD/waplkd_actd_variables']
    feature_all, file_dirs, dict_inputfile, dict_image, dict_sift, trainingset_filenames = util.read_features(
                    input_files, training_set=1)
    # util.save_features_small(input_files)

print '\n\n _______START FULL ANALYSIS_______ '; ininini_time = time.time()
print ' _________________________________'

if run['image_processing']:
    inini_time = time.time()
    for kk, file_dir in enumerate(file_dirs):  # each file_dir corresponds to one type of experimental data
        feature_all = []  # list of all features obtained for each training image in file_dirs
        if not is_dataset and (kk > 0): break
        print '\n\n\n _____________analysis experiment:', exp_names[kk], '_______'; iniexp_time = time.time()
        for jj, fileid in enumerate(sorted(os.listdir(file_dir)), start=0):
            if not is_dataset and (jj > 0): break
            # if jj > 1: break
            elif is_dataset: file_name = fileid.split(fileExt)[0]

            print '\n\n _______ file no.%d _______ ' % jj
            # # ============== INPUT PARAMETERS ========
            # # ========================================
            analysis_pixel_size = 10   # 10 [nm] <<< 160 [nm] (STORM res.) -> scale pixel size anal.p.s/ori.p.s
            dict_inputfile = dict(file_dir=file_dir, file_name=file_name, fileExt=fileExt, is_storm=is_storm,
                                                      out_channel='all',
                                  ispp=1, compute_ROI=0, crop=1, crop_range=[190, 250, 100, 170],
                                          pixelate=0,
                                          tessellate=1,
                                  original_pixel_size=160, photonconv=0.14, resolution=analysis_pixel_size)   # [nm]/[pixel], e.g. STORM
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
                    if not dict_inputfile['crop']:
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
                        print("\tDone roi (time =  %.2f seconds)" % (time.time() - start_time))
                    else:
                        crop_range = dict_inputfile.get('crop_range')
                        points_roi = util.points_2dcrop(data.points, crop_range)
                        # file_name_crop = file_name + '_'
                        # np.savetxt(output_dir + file_name_crop + '.txt', points_roi)
                        # points_roi = iproc.points_2dcrop(data.points2, crop_range)
                        file_name_crop = file_name + '_' + \
                                        str(crop_range[0])+str(crop_range[1]) + '_' +\
                                        str(crop_range[2]) + str(crop_range[3])
                        np.savetxt(output_dir + file_name_crop + '.txt', points_roi)
                        vproc.plot_points(points_roi)
                        # plt.savefig(output_dir + file_name + '_pp_roi' + '.pdf', bbox_inches='tight', dpi=400)
                        # percent = 80; points_reduced = util.remove_points_rnd(points, percent=percent)
                        # np.savetxt(output_dir + file_name + '_reducedpercent_' + str(percent) + '.txt', points_reduced)
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

                    print 'Computing Voronoi tessellation... '
                    vor = Voronoi(points_roi)  # compute Voronoi tessellation

                    print 'Computing Voronoi descriptors...'
                    vproc.compute_parameters(vor, dict_inputfile)   # new attribute in vor object: vor.areas
                    print "\tnuclear area = %.0f" % float(vor.areas_total), '[o. pixel size]2'

                    print 'Converting Voronoi tessellation into image (i.e., interpolate)...',
                    image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
                                                        interpolate_method=dict_image.get('interpolate_method'),
                                                        fill_value=0.0,
                                                        density_transform=dict_image['detect_densitytransform'])  # log 4 detect
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

                    print ("Done voronoi (time =  %.2f seconds)" % (time.time() - ini_time))

            # # ====== IMAGE PROCESSING ===================
            # # ===========================================
            print '\n _______IMAGE PROCESSING_______'; ini_time = time.time()
            plt.close('all')

            print 'Intensity-dependent computations...'
            dict_sift = dict(scale_pixel_size=scale_pixel_size, resolution=1,
                             original_pixel_size=dict_inputfile.get('original_pixel_size'),

                             # feature detection
                             t=10, feature_name='blob',  # 0.8, scale_ini=80
                             thresholding=1, threshold_percent=0.9, threshold_value=None, num_features=1000,
                             scale_range_is='nm',
                             scale_ini=20,
                             scale_end=1200,  # diam. of search, if pixel->analysispixel
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

            print 'Detecting and describing features (%s)...' % dict_sift.get('feature_name'); start_time = time.time()
            feature = iproc.find_feature(image, dict_sift, image_descr=image_descr)
            print ("\tDone (total time =  %.2f seconds)" % (time.time() - start_time))
            feature['nnd'] = stat.nnd_feature(feature, pixel_size=dict_image['scale_pixel_size']*dict_image['original_pixel_size'])
            if dict_inputfile.get('ispp'):
                feature['number_localizations'] = vproc.localizations_feature(vor, feature, dict_sift)
                feature['density'] = feature['number_localizations']/(np.pi*(0.5*feature['diameter'])**2)
                feature['clusterincluster'] = stat.count_clusterincluster(feature, analysis_pixel_size)
            print "\tnumber of (thresholded) features detected = %d" % feature.get('argmaxgrad')[0].shape

            print 'Plotting intensity-dependent Voronoi features...'
            start_time = time.time()
            iproc.plot_feature(image, feature, cmap='jet', norm='linear', plot_axis='on', blob_color='strength')
            plt.savefig(output_dir + file_name + '_features_image.pdf', bbox_inches='tight')
            print("Done (time =  %.2f seconds)" % (time.time() - start_time))

            # if dict_inputfile.get('ispp'):
            #     print 'Plotting intensity-dependent Voronoi features on Voronoi...',
            #     start_time = time.time()
            #     vproc.plot_feature(vor, feature, dict_sift, show_points=True, cmap='jet', norm='log', plot_axis='on')
            #     plt.savefig(output_dir + file_name + '_features_pp.pdf', bbox_inches='tight')
            #     print ("Done (time =  %.2f seconds)" % (time.time() - start_time))

            feature['file_name'] = file_name; feature['file_dir'] = file_dir
            # feature['image'] = image; feature['vor'] = vor
            feature_all.append(feature)

            # dict_output = dict(image_area=image.shape[0]*image.shape[1]*analysis_pixel_size**2,
            #                    ntotal=points.shape[0],
            #                    number_clusters=feature.get('argmaxgrad')[0].shape[0],
            #                    cluster_diameters=feature['diameter'].tolist(),
            #                    number_localizations_cluster=feature['number_localizations'].tolist(),
            #                    )
            # util.saveparameters(output_dir + file_name + '_output.txt', dict1=dict_output)
            print ("Done image processing (time =  %.2f seconds)" % (time.time() - ini_time))

            if not is_dataset:
                print 'Computing feature statistics of one image...'
                savefile_suffix = output_dir + feature['file_name']  # filename
                stat.statistic_descrip(feature, file_dirs=file_dirs[kk], ispp=dict_inputfile['ispp'],
                                       pixel_size=dict_image['scale_pixel_size']*dict_image['original_pixel_size'],
                                       savefile_suffix=savefile_suffix,
                                       radius=0, area=1, density=1, num_loc=1, nnd=1, density_cluster=1, area_voronoi=0,
                                       cluster_cluster=1)
        # if jj > 0:  # dir with more than one file analyzed
        #     savefile_suffix = output_dir + exp_names[kk]  # name current filedir
        #     print 'Computing feature statistics of dataset in %s',  savefile_suffix, '...'
        #     stat.statistic_descrip(feature_all, file_dirs=file_dirs[kk], ispp=dict_inputfile['ispp'],
        #                            pixel_size=analysis_pixel_size, savefile_suffix=savefile_suffix,
        #                            radius=0, area=1, density_cluster=1, area_voronoi=0, num_loc=1,
        #                            density=1, nnd=1, cluster_cluster=1)
        util.saveparameters(output_dir + exp_names[kk] + '_parameters.txt', dict1=dict_inputfile, dict2=dict_image,
                            dict3=dict_sift)
        file_variables = open(output_dir + exp_names[kk] + '_variables', 'wb')
        pickle.dump([feature_all, file_dirs[kk], dict_inputfile, dict_image, dict_sift], file_variables)
        file_variables.close()
        print ("Done experiment (time =  %.2f seconds)" % (time.time() - iniexp_time))
    if (kk > 0) and (kk < 2):  # 2 experiment
        feature_all = []
        for exp in exp_names:
            file_variables = open(output_dir + exp + '_variables', 'rb')
            [feature, file_dirs, dict_inputfile, dict_image, dict_sift] = pickle.load(file_variables)
            file_variables.close()
            feature_all.append(feature)
        savefile_suffix = output_dir + ''.join(exp_names)
        print 'Computing feature statistics of 2 datasets in', savefile_suffix, '...'
        stat.statistic_descrip(feature_all, file_dirs=file_dirs[0:2], ispp=dict_inputfile['ispp'],
                               pixel_size=dict_image['scale_pixel_size']*dict_image['original_pixel_size'],
                               savefile_suffix=savefile_suffix,
                               radius=0, area=1, density=1, num_loc=1, nnd=1, density_cluster=1, area_voronoi=1,
                               cluster_cluster=1, strength=1)
        # same graph as before with boxplots
        labels = np.asarray([np.where(feature['file_dir'] == np.asarray(file_dirs))[0][0] for feature in feature_all
                                  for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                  enumerate(orientation)])  # trick
        stat.words_statistic_descrip(feature_all, labels, cmap=mpl.colors.ListedColormap(['black', 'red']),
                                     xlabel=r'circle $\quad$rectangle', #r'hFb $\quad$TSA-hFb',
                                     pixel_size=dict_sift['original_pixel_size'] * dict_sift['scale_pixel_size'],
                                     savefile_suffix=savefile_suffix, pt='_experiments',
                                     radius=0, area=1, num_loc=1, density=1, nnd=1, strength=1, experiment=1,
                                     cluster_density=1, cluster_cluster=1)
        # # boxplots for one dir and multiple images
        # stat.words_statistic_descrip(feature_all, np.asarray([ff for ff, feature in enumerate(feature_all)
        #                              for ii, orientation in enumerate(feature['orientation']) for jj, ori in
        #                              enumerate(orientation)]), cmap=mpl.colors.ListedColormap(['black', 'black']),
        #                              xlabel='',
        #                              pixel_size=dict_sift['original_pixel_size'] * dict_sift['scale_pixel_size'],
        #                              savefile_suffix=savefile_suffix,
        #                              radius=1, area=1, num_loc=1, density=1, nnd=1, strength=1, experiment=0,
        #                              cluster_density=1, cluster_cluster=1)
        util.saveparameters(savefile_suffix + '_parameters.txt', dict1=dict_inputfile, dict2=dict_image, dict3=dict_sift)
        file_variables = open(savefile_suffix + '_variables', 'wb')
        pickle.dump([feature_all, file_dirs, dict_inputfile, dict_image, dict_sift], file_variables); file_variables.close()
    print ("Done (total time =  %.2f seconds)" % (time.time() - inini_time))

# # ====== IMAGE ANALYSIS ====================
# # ==========================================
feature_all_original = feature_all  # save original features before filtering etc.
# feature_all = feature_all_original
dict_analysis = dict(filter=0, name=['blob_size'],
                               ball_radius_percentile=5, diameter_range=[300, 1000], word_id=[1], min_numloc=3,
                               experiment_name=file_dirs[0], threshold_percent=0.7, threshold_value=None,
                     reduce_dimension=1, method_rd='class_equalsize')
if dict_analysis['filter']:
    try: siftclusters
    except NameError: siftclusters = None
    feature_all_filtered = feature_all; siftclusters_filtered = siftclusters
    if (dict_analysis['name'] is 'word_id') or (dict_analysis['name'] is 'histogram_distance') and siftclusters is None:
        print 'error: run image_analysis before filtering by labels (word ids)'
    for ii, fn in enumerate(dict_analysis['name']):  # filter detected features by different criteria
        print 'filtering according to', fn
        feature_all_filtered, feature_pos_filtered, siftclusters_filtered = \
                        stat.filter_features(fn, feature_all_filtered, dict_analysis, siftclusters=siftclusters_filtered)
        if ii == 0: feature_pos_filtered_prev = feature_pos_filtered; continue
        feature_pos_filtered_prev = feature_pos_filtered_prev[feature_pos_filtered]
    feature_pos_filtered = feature_pos_filtered_prev
    feature_all = feature_all_filtered

if run['image_analysis']:
    plt.close('all')
    reload(stat); reload(iproc); reload(vproc); reload(util)
    # # # ==== CREATE VOCABULARY, unsupervised =====
    if dict_analysis['reduce_dimension']:
        feature_all_reduced, feature_pos_reduced = stat.reduce_dimensionality(feature_all, method=dict_analysis[
            'method_rd'], n_cluster=int(0.5*iproc.number_features(feature_all)))
        feature_all_filtered = feature_all_reduced; feature_pos_filtered = feature_pos_reduced
        feature_all = feature_all_reduced
    if dict_sift.get('compute_orientation'):
        print '\n _______ IMAGE ANALYSIS_______'; inini_time = time.time()
        k0 = 4; kn = 4; sserror = []; avsilhouette = []; percentilesilhouette = []
        for k in range(k0, kn+1, 1):
            print 'Creating vocabulary with %d words...' % k
            siftclusters, histogram_weighted = stat.create_vocabulary('kmeans', feature_all, dict_sift,
                                                                      n_cluster=k, weight=0)
            # # # #  plot the feature class with k+1 according to the colorcode at k. Comment if it feels unnecessary.
            # # # if k == k0: centers_permuted0 = siftclusters.cluster_centers_
            # # # if k > k0: centers_permuted0, siftclusters = util.permute_labels(siftclusters, k, centers_permuted0)
            #
            if len(file_dirs) > 1:
                savefile_suffix_all = output_dir + ''.join(exp_names) + '_siftclusters%d' % k
            else:
                savefile_suffix_all = output_dir + feature_image['file_name'] + '_siftclusters%d' % k
            stat.scatterplot_vocabulary(feature_all, siftclusters, n_cluster=k, cmap=util.discrete_cmap(k, "jet"),
                                        savefile_suffix=savefile_suffix_all, cmap_charact=None,
                                        file_dirs=file_dirs,
                                        pixel_size=dict_sift['original_pixel_size'] * dict_sift['scale_pixel_size'],
                                        histogram=histogram_weighted)
            # # # # # stat.sift2shape(siftclusters, cmap=util.discrete_cmap(k, "jet"), savefile_suffix=savefile_suffix_all)
            stat.words_statistic_descrip(feature_all, siftclusters.labels_, cmap=util.discrete_cmap(k, "jet"),
                                         pixel_size=dict_sift['original_pixel_size']*dict_sift['scale_pixel_size'],
                                         ylabel=exp_names, file_dirs=file_dirs,
                                         savefile_suffix=savefile_suffix_all,
                                         radius=0, area=1, num_loc=1, density=1, nnd=0, strength=1,
                                         experiment=1, cluster_density=0, cluster_cluster=0)
            # stat.codewords_hist(feature_all, siftclusters)  # image representation
            # # # # plots and statistics on one selected image, optional (change siftclusters.labels_ ->
            # # siftclusters_labels_image):
            # image_no = 30; feature_image = feature_all[image_no]
            # # image=feature_image['image']; vor=feature_image['vor']
            # dict_inputfile['file_dir'] = feature_image['file_dir']
            # dict_inputfile['file_name'] = feature_image['file_name']
            # savefile_suffix = output_dir + feature_image['file_name'] + '_siftclusters%d' % k  # or None
            # if dict_inputfile.get('ispp'):
            #     data = util.pointpattern(); data.read(dict_inputfile, save_output_dir=output_dir, plot=0)
            #     points = data.points; vor = Voronoi(points)
            # siftclusters_labels_image = util.select_labels_image(feature_all, siftclusters.labels_, image_no=image_no)
            # # vproc.compute_parameters(vor, dict_inputfile)  # new attribute in vor object: vor.areas
            # # image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
            # #                                     interpolate_method=dict_image.get('interpolate_method'),
            # #                                     fill_value=0.0,
            # #                                     density_transform=dict_image['detect_densitytransform'])
            # # iproc.plot_feature(image, feature_image, cmap='jet', norm='linear', plot_axis='on',
            # #                    blob_color='class', ori_color=siftclusters_labels_image,
            # #                    ori_cmap=util.discrete_cmap(k, "jet"))
            # # plt.savefig(savefile_suffix + '_features_image.pdf', bbox_inches='tight')
            # if dict_inputfile['ispp']:
            #     vproc.plot_feature(vor, feature_image, dict_sift, norm='log', blob_color='class',
            #                        ori_color=siftclusters_labels_image, ori_cmap=util.discrete_cmap(k, "jet"))
            #     plt.savefig(savefile_suffix + '_features_pp.pdf', bbox_inches='tight')

            # # # # filtered PC points:
            # # feature_all_filtered, feature_pos_filtered, siftclusters_filtered = \
            # #    stat.filter_features('strength', feature_all, siftclusters=siftclusters,
            # #                        ball_radius_percentile=2, word_id=1,
            # #                        diameter_range=[100, 1000], experiment_name=file_dirs[0], threshold_percent=0.1)
            stat.scatterplot_vocabulary(feature_all, siftclusters, n_cluster=k, cmap=util.discrete_cmap(k, "jet"),
                                        savefile_suffix=savefile_suffix_all+'_filtered',
                                        filter_pos=feature_pos_filtered)
            # stat.words_statistic_descrip(feature_all_filtered, siftclusters_filtered.labels_,
            #                              cmap=util.discrete_cmap(k, "jet"), ylabel=r'ActD$\quad$DMSO',
            #                              pixel_size=dict_sift['original_pixel_size'] * dict_sift['scale_pixel_size'],
            #                              savefile_suffix=savefile_suffix_all+'_filtered',
            #                              radius=1, area=1, num_loc=1, density=1, nnd=1, strength=1,
            #                              experiment=1, cluster_density=1, cluster_cluster=1)
            # siftclusters_labels_image = util.select_labels_image(feature_all_filtered, siftclusters_filtered.labels_, image_no=image_no)
            # iproc.plot_feature(image, feature_all_filtered[image_no], cmap='jet', norm='linear', plot_axis='on',
            #                    blob_color='class', ori_color=siftclusters_labels_image, ori_cmap=util.discrete_cmap(k, "jet"))
            # plt.savefig(savefile_suffix + '_features_image_filtered.pdf', bbox_inches='tight')
            if dict_inputfile['ispp']:
                vproc.plot_feature(vor, feature_all_filtered[image_no], dict_sift, norm='log', blob_color='class',
                                   ori_color=siftclusters_labels_image, ori_cmap=util.discrete_cmap(k, "jet"))
                plt.savefig(savefile_suffix + '_features_pp_filtered.pdf', bbox_inches='tight')

            if hasattr(siftclusters, 'inertia_'): sserror.append(siftclusters.inertia_)
            silhouette, silhouette_percentile = stat.compute_silhouette(histogram_weighted, siftclusters.labels_,
                                         savefile_suffix=savefile_suffix_all, cmap=util.discrete_cmap(k, "jet"))
            avsilhouette.append(silhouette); percentilesilhouette.append(silhouette_percentile)

        fig2, ax2 = plt.subplots()
        ax2.plot(range(k0, kn+1), sserror, 'k.-', markersize=10)
        ax2.set_ylabel('total within sum of squares'); ax2.set_xlabel(r'number of clusters $k$')
        fig2.savefig(output_dir + ''.join(exp_names) + '_withinsumofsquares.pdf', bbox_inches='tight')

        # fig3, ax3 = plt.subplots()
        # ax3.plot(range(k0, kn+1), avsilhouette, 'k.-', markersize=10)
        # ax3.set_ylabel(r'average silhouette coefficient $\overline{s(i)}$'); ax3.set_xlabel(r'number of clusters $k$')
        # fig3.savefig(output_dir + ''.join(exp_names) + '_avsilhouette.pdf', bbox_inches='tight')
        #
        # fig4, ax4 = plt.subplots()
        # ax4.plot(range(k0, kn+1), 100*percentilesilhouette, 'k.-', markersize=10)
        # ax4.set_ylabel(r'clusters below $\overline{s(i)}$ [%]'); ax4.set_xlabel(r'number of clusters $k$')
        # fig4.savefig(output_dir + ''.join(exp_names) + '_percentilesilhouette.pdf', bbox_inches='tight')
        # print ("Done (total time =  %.2f seconds)" % (time.time() - inini_time))

print ("\n\n***FINISHED ANALYSIS (total time =  %.2f seconds)" % (time.time() - ininini_time))

# # if __name__ == '__main__':
# #     main()
