import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams["text.usetex"] = True; matplotlib.rcParams['font.family'] = 'serif'  # configure latex plots
matplotlib.rcParams.update({'font.size': 16})


def plot_hist(data, bins=None, scale='lin', xlabel={}, cmap=None, num_bins=100, xticks=None, xtickslabel=None,
              alpha=1):
    """
    This function plots histograms. Range: min-max of data.

    Input:
    ----------------
    scale (string): log or lin xaxis of the histogram. Default = 'log'
            i.e., bins=np.logspace(np.log10(ini), np.log10(np.max(data)), num=num_bins)
                  bins=np.linspace(ini, np.max(data), num=num_bins)
    num_bins (float): number of bins. Default = 100
    bins (1d array): if it is not None, then bins=bins and 'lin' scale

    """
    fig, ax = plt.subplots()

    if scale is 'log':
        data = data[np.where((data != float('inf')) & (data > 0))]  # eliminate areas that are set to inf, -1 or 0
    else:
        data = data[np.where((data != float('inf')) & (data >= 0))]
    if bins is not None:
        n, bins, patches = ax.hist(data, bins=bins, histtype='bar',
                                    weights=np.zeros_like(data) + 1. / data.size, color='w')  # , color='k')
        if cmap is not None:
            for c, p in zip(range(len(bins)-1), patches):
                plt.setp(p, 'facecolor', cmap(c), alpha=alpha)
        if xticks is not None:
            if xtickslabel is None: xtickslabel=xticks
            plt.xticks(xticks, xtickslabel)
        fig.show()
    else:
        ini = np.min(data)
        if scale is 'log':
            ax.hist(data, bins=np.logspace(np.log10(ini), np.log10(np.max(data)), num=num_bins),
                                   histtype='step', weights=np.zeros_like(data) + 1. / data.size,
                                   color='k', alpha=alpha)
            ax.set_xscale("log")
        else:
            ax.hist(data, bins=np.linspace(ini, np.max(data), num=num_bins), histtype='step',
                                   weights=np.zeros_like(data) + 1. / data.size, color='k', alpha=alpha)
            # ax.hist(data, bins=num_bins, histtype='step',
            #                        weights=np.zeros_like(data) + 1. / data.size, color='k')
    plt.ylabel(r'frequency'); plt.xlabel(xlabel); ax.hold(0)
    ax.set_ylim(0, 0.6)
    # ax.hist(densities, bins='rice', histtype='step',  color='k'); plt.ylabel(r'counts')

    return fig, ax


def plot_boxplot(data, scale='lin', bptype='violin', xticklabels='', xlabel='', ylabel='values', cmap=None,
                 widths=0.95, alpha=0.5):
    """
    This function plots boxplot.

    Input:
    ----------------
    data: list of 1 or more arrays (depending on the number of boxplots per figure), e.g. [np.array()] or
            [np.array(), np.array(),...]
    xlabel: list of strings corresponding to the 1 or more boxplots in the figure, e.g. xlabel=['control'] or [
            'control', 'experiment']

    Output:
    --------------
    result (dict): cmeans, cmedians, ... see: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.violinplot.html
    """

    if not isinstance(data, list): data = [data]  # correct if not list input data

    if scale is 'log':
        data = [bp[np.where((bp != float('inf')) & (bp > 0))] for bp in data] # eliminate areas that are set to inf, -1 or 0
    else:
        data = [bp[np.where((bp != float('inf')) & (bp >= 0))] for bp in data]

    if scale is 'log':
        values = [np.log10(data)]
        ylabel = ylabel + '$\quad$' + r'[$\log_{10}$]'
    else: values = data

    fig, ax = plt.subplots()
    if bptype == 'violin':   # plot violin plot
        result = ax.violinplot(values, widths=widths, showmeans=1, showmedians=True, showextrema=True)
        for ii, pc in enumerate(result['bodies']):
            pc.set_alpha(alpha)
            if cmap is None:
                pc.set_facecolor('gray')
                pc.set_linewidth(1)
            else:
                pc.set_facecolor(cmap(ii))
                pc.set_linewidth(1)
        for pc in [result[ii] for ii in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans')]:
            pc.set_edgecolor('black')
            pc.set_linewidth(2)
            pc.set_alpha(alpha)
            if pc == result['cmeans']:
                pc.set_linestyle('--')
    else:  # plot box plot
        result = ax.boxplot(values)

    # # adding horizontal grid lines
    # ax.yaxis.grid(True)
    # ax.set_xticks([y + 1 for y in range(len(values))])
    # ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)

    if xticklabels == '': xticklabels=[str(i+1) for i in range(len(values))]
    plt.setp(ax, xticks=[y + 1 for y in range(len(values))], xticklabels=xticklabels, ylabel=ylabel, xlabel=xlabel)
    plt.show()

    return result


def sample_statistics(data):

    mean = round(np.nanmean(data), 2)
    std = round(np.nanstd(data, ddof=1), 2)  # sum()/(N-ddof)
    median = round(np.nanmedian(data), 2)

    return {'mean': mean, 'std': std, 'median': median}


def statistic_descrip(feature_all, ispp=1, pixel_size=1, savefile_suffix=None):
    """
    plot and save graphs regarding detection of clusters (diameter, densities...)

    Input:
    ----------------
    feature_all: it can be a list of features extracted from a data set or just one feature from one sample.
    """

    if not isinstance(feature_all, list): feature_all = [feature_all]  # only one sample

    blob_diameters_all = np.asarray([diameter for feature in feature_all for diameter in feature['diameter']])
    cluster_density_all = np.asarray([feature['argmaxgrad'][0].shape[0] / float(pixel_size**2 *
                                      feature['image'].shape[0] * feature['image'].shape[1]) for feature in
                                      feature_all])
    # nnd_localizations = iproc.nnd_feature(feature, dict_sift)
    plot_hist(0.5*blob_diameters_all, scale='log', num_bins=np.unique(feature.get('scale')).shape[0]+1,
              xlabel=r'blob radius R [nm]')
    plt.savefig(savefile_suffix + '_blobradius_hist.pdf', bbox_inches='tight')
    plot_boxplot(0.5 * blob_diameters_all, bptype='violin', ylabel=r'blob radius R [nm]')
    plt.savefig(savefile_suffix + '_blobradius_boxplot.pdf', bbox_inches='tight')
    # plot_hist(np.pi*(blob_diameters_all/2.)**2, num_bins=50, xlabel=r'blob area [nm2]')
    # plot_boxplot(np.pi * (blob_diameters_all / 2.) ** 2, bptype='violin', ylabel=r'blob area [nm2]')
    # plt.savefig(savefile_suffix + '_blobareas_boxplot.pdf', bbox_inches='tight')
    # errorbar_featureresponse(feature, dict_sift, xlabel=r'blob diameter [nm]')
    # print '\tNumber of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0]
    # print '\tDensity of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0] / float(pixel_size *
    #                 feature_all[0]['image'].shape[0] * feature_all[0]['image'].shape[1]), '[cluster/nm2]'
    plot_boxplot(cluster_density_all, scale='lin', bptype='violin', ylabel=r'$\kappa$ [clusters/nm$^2$]')
    plt.savefig(savefile_suffix + '_blobdensitieskappa_boxplot.pdf', bbox_inches='tight')
    if ispp:
        density_all = np.asarray([density for feature in feature_all for density in feature['density']])
        vor_areas_all = np.asarray([area for feature in feature_all for area in feature['vor'].areas])
        num_loc_all = np.asarray([num_loc for feature in feature_all for num_loc in feature['number_localizations']])
        # percentage_number_localizations = np.asarray([numloc / points.shape[0] * 100 for feature in feature_all for
                                            #  numloc in feature['number_localizations']])
        # vareas_statistics = sample_statistics(vor.areas[vor.areas < float('inf')]*
        #                                       dict_inputfile['original_pixel_size']**2)
        # densities_statistics = sample_statistics(density_all)
        # feature_statistics = sample_statistics(feature['number_localizations'])
        # print '\tPercentage of localizations in clusters:\t', np.sum(percentage_number_localizations), '%'
        # print '\tVoronoi polygon areas: \t', vareas_statistics, '[nm2]'
        # print '\tNumber of loc. per blob: \t', feature_statistics, '[loc/blob]'
        # print '\tLocal blob densities: \t', densities_statistics, ' [loc/nm2]'
        # plot_hist(vor_areas_all, scale='log', num_bins=50, xlabel=r'Voronoi polygon area [nm$^2$]')
        plot_boxplot(vor_areas_all, scale='log', bptype='violin', ylabel=r'Voronoi polygon area [nm$^2$]')
        plt.savefig(savefile_suffix + '_voronoiareas_boxplot.pdf', bbox_inches='tight')
        plot_boxplot(num_loc_all, scale='log', bptype='violin', ylabel=r'N$^{cluster}$ [points/cluster]')
        plt.savefig(savefile_suffix + '_blobnumloc_boxplot.pdf', bbox_inches='tight')
        plot_hist(num_loc_all, scale='log', num_bins=50, xlabel=r'N$^{cluster}$ [points/cluster]')
        plt.savefig(savefile_suffix + '_blobnumloc_hist.pdf', bbox_inches='tight')
        # plot_hist(densities, scale='lin', num_bins=50, xlabel=r'blob density [localizations/nm$^2$]')
        plot_boxplot(density_all, bptype='violin', ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
        plt.savefig(savefile_suffix + '_blobdensities_boxplot.pdf', bbox_inches='tight')


def errorbar_featureresponse(feature, dict_sift, xlabel={}):
    """
    Plot the error bar strength response vs blob diameter

    Input:
    -------------------
    feature (dictionary)
    dict_analysis (dictionary)

    """
    strength = feature.get('strength')
    scale = feature.get('scale')
    scale_range = feature.get('scale_range')
    analysis_pixel_size = dict_sift.get('scale_pixel_size')*dict_sift.get('original_pixel_size')

    err = np.zeros(shape=scale_range.shape)
    y = np.zeros(shape=scale_range.shape)

    for ii in range(len(scale_range)):
        pos_equal_scale = np.where(scale == scale_range[ii])
        str = strength[pos_equal_scale]
        y[ii] = np.median(str)
        if math.isnan(y[ii]): y[ii] = 0
        err[ii] = np.std(str)

    x = 3*np.sqrt(scale_range)*analysis_pixel_size
    plt.figure(); plt.errorbar(x, y, err, fmt='ko')
    plt.xlabel(xlabel); plt.ylabel('max$_t\{\Delta_{\gamma-norm}\}$')


def compute_rms_deviation(points, area, width, aspect_ratio, bg, kwargs, plot=False):
    """
    Deviations of localization estimations. Based on (Mortensen et al., 2009), 'Optimized localization analysis for
    single-molecule tracking and super-resolution microscopy', Eq.(6). Sigma is defined in (Thompson et al., 2002).
    Computes the theroetical r.m.s. deviations of estimatiors = the square root of twice th expressions for the variance
    of our estimates, mu_x. See also (Insight3 manual Joe, page 3).
    """
    a = kwargs.get('original_pixel_size', 160)
    pc = kwargs.get('photonconv', 0.14)  # nm/pixel

    b = bg*pc
    s = 0.5*width*np.sqrt(aspect_ratio)
    photons = pc*area  # number of photons for each molecule in the molecule list. The photon conversion factor,
    #  eg. 0.41 for STORM and 0.14 for NSTORM
    N = photons  # total photon count, number of photons in the image
    mu_x = (s**2 + a**2/12.)/N*(16./9 + 8.*np.pi*s**2*b**2/(N*a**2))
    rms = np.sqrt(2*mu_x)

    if plot is True:
        import matplotlib.colors as colors
        import matplotlib.cm as cmx

        cmap = plt.cm.gray
        cNorm = colors.Normalize(vmin=0, vmax=np.max(rms))  # vmax=values[-1]) . LogNorm, Normalize
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)  # norm=cNorm
        scalarMap.set_array(rms)

        for ii in range(len(points)):
            deviation = rms[ii]
            blob_color = scalarMap.to_rgba(deviation)
            ax = plt.plot(points[ii][0], points[ii][1], color=blob_color, markersize=1)
            plt.hold(True)
        plt.colorbar(scalarMap, label='$\sqrt{2\mu_x}$')
        plt.hold(False)

    return rms


def scatterplot_vocabulary(feature_all, kmeans, n_cluster=3, cmap=None, savefile_suffix=None):
    """
    This functions plots in the PC1-PC2 space the clustered distribution of all detected features to create the
    vocabulary.
    """
    from sklearn.decomposition import PCA

    if not isinstance(feature_all, list): feature_all = [feature_all]

    histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
    histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
    if len(histogram) == 0:
        print 'error: please, compute orientation and descriptors histograms for SIFT'
        return None

    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_  # kmeans.cluster_centers_[0] = centroid of cluster 0

    print '\tcomputing PCA for visualization...'
    pca = PCA(n_components=2)
    X = np.append(histogram, cluster_centers, axis=0)  # shape (n_samples, n_features)
    X_reduced = pca.fit_transform(X)  # array-like, shape (n_samples, n_components)
    # histogram_reduced = pca.fit_transform(histogram)  # array-like, shape (n_samples, n_components)
    # cluster_centers_reduced = pca.fit_transform(cluster_centers)  # array-like, shape (n_samples, n_components)
    # print cluster_centers_reduced

    # fig, ax = plt.subplots()
    plt.figure()
    if cmap is None: cmap = plt.cm.get_cmap('Set1')
    plt.scatter(X_reduced[0:histogram.shape[0], 0], X_reduced[0:histogram.shape[0], 1], c=labels,
                alpha=0.7, cmap=cmap)
    plt.hold(True)
    # aux = 0
    # orientation = feature.get('orientation')
    # argmaxgrad = feature.get('argmaxgrad')
    # for i in range(argmaxgrad[0].size):
    #     for ii, ori in enumerate(np.asarray(orientation[i])):  # if empty => no in X_reduced
    #         # ax.annotate('(%d,%d), ori=%.1f' % (argmaxgrad[0][i], argmaxgrad[1][i], ori),
    #         #             xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
    #         ax.annotate('%d' % labels[aux], xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
    #         aux += 1
    plt.scatter(X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 0],
               X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 1], marker='*', s=10**2,
               color=[cmap(ii) for ii in range(n_cluster)])
    # plt.plot(cluster_centers_reduced[:, 0], cluster_centers_reduced[:, 1], 'k*', markersize=10)
    # for ii in range(n_cluster):
    #     ax.annotate('%d' % ii, xy=(X_reduced[histogram.shape[0] + ii, 0],
    #                                X_reduced[histogram.shape[0] + ii, 1]), size=40)
    # ax.annotate('%.d', (X_reduced[hist_ind, s0], X_reduced[hist_ind, 1]), size=10)
    plt.xlabel(r'PC$_1$'); plt.ylabel(r'PC$_2$')  # ;plt.title('k=%d' % n_cluster)
    plt.show(); plt.hold(False)

    if savefile_suffix is not None: plt.savefig(savefile_suffix + '.pdf', bbox_inches='tight')

    # arg = 0; hist_ind = 0
    # for hist in histogram_descr:
    #     if not hist:  # skip that blob, out of range
    #         arg += 1
    #         continue
    #     bx = argmaxgrad[0][arg]; by = argmaxgrad[1][arg]
    #     t = scale[arg]
    #     arg += 1
    #     for jj in range(len(hist)):
    #         ax.annotate('%.0f,%.0f' % (bx, by), (X_reduced[hist_ind, 0]-0.01, X_reduced[hist_ind, 1]-0.01), size=10)
    #         ax.annotate('%.0f' % t, (X_reduced[hist_ind, 0], X_reduced[hist_ind, 1]), size=10)
    #         hist_ind += 1


def create_vocabulary(feature_all, kwargs_sift={}, n_cluster=3, init="k-means++"):
    """
    This function

    Input:
    ---------
    feature_all: list of feature (argmaxgrad, orientation, scale, histogram_descr, ...) of all training images

    Output:
    ---------
    kmeans: attributs are labels_, cluster_centers_

    """
    from sklearn.cluster import KMeans

    if not isinstance(feature_all, list): feature_all = [feature_all]

    histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
    histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])  # some hist can be []
    # print 'size all histograms: ', len(histogram)
    if len(histogram) == 0:
        print 'error: please, compute orientation and descriptors histograms for SIFT'
        return None

    kmeans = KMeans(n_clusters=n_cluster, random_state=0, init=init).fit(histogram)

    # labels = kmeans.labels_
    # cluster_centers = kmeans.cluster_centers_  # kmeans.cluster_centers_[0] = centroid of cluster 0
    # n_bins_descr = kwargs_sift.get('n_bins_descr', 4)
    # n_hist = kwargs_sift.get('n_hist', 4)

    # # plot
    # for jj in range(cluster_centers.shape[0]):
    #     plt.figure()
    #     max_ori = np.max(cluster_centers[jj])
    #     for ii in range(0, n_hist * n_hist):
    #         ind = 4 * (3 - ii // 4) + ii % 4
    #         plt.subplot(4, 4, ii + 1)
    #         plt.bar(np.arange(n_bins_descr), cluster_centers[jj][n_bins_descr * ind:n_bins_descr * ind +n_bins_descr],
    #                 color='0.5',
    #                 align='center')
    #         plt.xticks(np.linspace(0, n_bins_descr, 5), ['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
    #         plt.ylim([0, max_ori])
    #         plt.hold(True)

    return kmeans


def words_characteristics(feature_all, labels, cmap, pixel_size=1, savefile_suffix=None,
                          blob_diameter_thr=[-float('inf'), float('inf')]):
    """
    This function plots histograms and violin plots of the scales, densities, etc., corresponding to each cluster in the
    vocabulary

    Input:
    ---------------------
    feature_all (list): vector of all feature computed from training data
    labels: from kmeans.labels_, output of the clustering of the trained features, each point corresponds to the
            orientation of a features of a given training image key() = labels_, cluster_centers_
    blob_diameters_thr [a,b]: results for the detected blobs of a particular scale range (diameter)

    Alternatively, features_all[n] & kmeans_labels_image = kmeans.labels_[~n]

    Output:
    ---------------------

    """
    if not isinstance(feature_all, list): feature_all = [feature_all]

    # words histograms
    labels_sub = []; hist_ind = 0
    for feature in feature_all:
        blob_diameters = 3*np.sqrt(feature['scale'])*pixel_size
        orientation = feature.get('orientation', [])
        if len(orientation) == 0:
            print 'Error: SIFT descriptors missing.'
            return -1
        for ii, ori in enumerate(orientation):
            for jj in ori:  # loop just to update hist_ind
                if blob_diameter_thr[0] < blob_diameters[ii] < blob_diameter_thr[1]:
                    labels_sub.append(labels[hist_ind])  # if no thr & no limits = labels
                hist_ind += 1
    plot_hist(np.asarray(labels_sub),
              bins=np.append(np.unique(labels_sub), np.max(labels_sub)+1)-0.5,
              cmap=cmap, xlabel=r'visual word - cluster num.', xticks=np.unique(labels_sub),
              xtickslabel=[str(e) for e in np.unique(labels_sub)+1], alpha=0.7)
    if savefile_suffix is not None:
        plt.savefig(savefile_suffix + '_BoWhistogram.pdf', bbox_inches='tight')

    # scale vs words ; densities vs words ; strength-laplacian vs words
    scale_all = [scale for feature in feature_all for scale in feature['scale']]
    number_localizations_all = [number_localizations for feature in feature_all for number_localizations in
                                 feature['number_localizations']]
    strength_all = [strength for feature in feature_all for strength in feature['strength']]
    scale_ori = np.asarray([scale_all[ii] for feature in feature_all
                            for ii, hist in enumerate(feature['histogram_descr']) for jj in range(len(hist))])
    number_localizations_ori = np.asarray([number_localizations_all[ii] for feature in feature_all
                                           for ii, hist in enumerate(feature['histogram_descr'])
                                           for jj in range(len(hist))])
    strength_ori = np.asarray([strength_all[ii] for feature in feature_all
                            for ii, hist in enumerate(feature['histogram_descr']) for jj in range(len(hist))])
    radius = []; densities = []; strength = []
    for l in np.unique(labels):
        pos = np.where(labels == l); rad = pixel_size*1.5*np.sqrt(scale_ori[pos])
        radius.append(rad)
        densities.append(number_localizations_ori[pos] / (rad ** 2 * np.pi))
        strength.append(strength_ori[pos])
    # plot boxplots
    plot_boxplot(radius, bptype='violin', xlabel=r'visual word - cluster num.', ylabel=r'blob radius R [nm]', cmap=cmap,
                 widths=0.8, alpha=0.7)
    if savefile_suffix is not None:
        plt.savefig(savefile_suffix + '_scaleVSwords.pdf', bbox_inches='tight')
    plot_boxplot(densities, bptype='violin', xlabel=r'visual word - cluster num.',
                 ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]', cmap=cmap, widths=0.8, alpha=0.7)
    if savefile_suffix is not None:
        plt.savefig(savefile_suffix + '_denstitiesVSwords.pdf', bbox_inches='tight')
    plot_boxplot(strength, bptype='violin', xlabel=r'visual word - cluster num.',
                 ylabel=r'strength - laplacian, max$_t\{\Delta_{\gamma-norm}\}$', cmap=cmap, widths=0.8, alpha=0.7)
    if savefile_suffix is not None:
        plt.savefig(savefile_suffix + '_strengthVSwords.pdf', bbox_inches='tight')

    # for ii, feature in enumerate(feature_all):
    #     # pos_empty = [ii for ii, hist in enumerate(feature_all[ii]['histogram_descr']) if hist == []]
    #     for label in np.unique(kmeans.labels_):
    #         scale_ori = [feature_all[ii]['scale'][ii] for jj in range(len(hist))
    #                      for ii, hist in enumerate(feature_all[ii]['histogram_descr'])]
    #         diameter = 3*np.sqrt(scale_ori[np.where(scale_ori == label)])


def sift2shape(kmeans, n_hist=16, n_bins_descr=8, cmap=None, savefile_suffix=None, alpha=0.7):
    """
    After unsupervised-clustered algorithm for vocabulary contruction, build histograms and/or shapes of the estimated
    words (or centroids)

    Input:
    -----------------

    """

    # matplotlib.rcParams.update({'font.size': 10})

    if cmap is None: cmap = plt.cm.get_cmap('Set1')
    # matplotlib.rcParams.update({'font.size': 10})

    for no, hist in enumerate(kmeans.cluster_centers_):
        max_ori = np.max(hist)
        fig, axes = plt.subplots(nrows=4, ncols=4, sharex='col', sharey='row'); # fig.suptitle(r'word %d' % (no+1))
        plt.setp(axes, xticks=np.linspace(0, n_bins_descr, 5), xticklabels=['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
        # plt.xticks(np.linspace(0, n_bins_descr, 5), ('0', '\pi/2', '\pi', '3\pi/2', '2\pi'))
        for ii in range(0, n_hist):
            ind = 4 * (3 - ii // 4) + ii % 4
            plt.subplot(4, 4, ii + 1)
            plt.bar(np.arange(n_bins_descr), hist[n_bins_descr * ind:n_bins_descr * ind + n_bins_descr],
                    align='center', alpha=alpha, color=cmap(no))  # color='0.5'
            plt.xticks(np.linspace(0, n_bins_descr, 5), ['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
            plt.tick_params(labelsize=10)
            plt.ylim([0, max_ori]); plt.hold(True)
        # fig.tight_layout()
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_sift2shapeWORD%d' % no + '.pdf', bbox_inches='tight')


def filter_features(feature_all, kmeans, ball_radius):
    """
    ...
    """
    import copy

    if not isinstance(feature_all, list): feature_all = [feature_all]

    histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
    histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr
                            for hist_sub in hist])  # some hist can be []
    if len(histogram) == 0:
        print 'error: please, compute orientation and descriptors histograms for SIFT'
        return None

    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_  # kmeans.cluster_centers_[0] = centroid of cluster 0

    # for ii, center in enumerate(cluster_centers):
    #     cluster_id_pos = np.where(labels == ii)
    #     cluster_in = np.where(np.linalg.norm(histogram[cluster_id_pos] - cluster_centers[0], axis=1) > ball_radius)
    #     feature_pos = cluster_id_pos[0][cluster_in[0]]
    cluster_id_pos = np.where(labels == 0)
    cluster_in = np.where(np.linalg.norm(histogram[cluster_id_pos] - cluster_centers[0], axis=1) > ball_radius)
    feature_pos = cluster_id_pos[0][cluster_in[0]]

    argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
    for ii, feature in enumerate(feature_all):
        for jj, ori in enumerate(feature['orientation']):
            for kk, o in enumerate(ori):  # only count if !=[]
                image_nos.append(ii)  # for each label/feature/ori position, what is the corresponding image
                argmaxgrad_feature_pos.append(jj)
                # for each feature position, position of argmaxgrad in that image
                ori_feature_pos.append(kk)  # for each feature position, position of ori in that orientation

    feature_all_filtered = copy.deepcopy(feature_all)
    # print feature_pos
    # print 'image_nos=', image_nos
    # print 'argmaxgrad_feature_pos=', argmaxgrad_feature_pos
    # print 'ori_feature_pos=', ori_feature_pos
    aux = 0; blob_loc_prev = None; blob_image_prev = None
    for pos in feature_pos:
        pos_ori = pos
        # print 'ori_feature_pos[pos_ori]=', ori_feature_pos[pos_ori], '; argmaxgrad_feature_pos[pos]=', argmaxgrad_feature_pos[pos]
        if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # same blob
            # print 'in'
            aux += 1; pos_ori -= aux
        else: aux = 0
        # print 'ori_feature_pos[pos_ori]=', ori_feature_pos[pos_ori]
        # print feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]]
        feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
            np.delete(feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]], ori_feature_pos[
            pos_ori])
        del feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
            pos_ori]]
        # print 'after deletion, ', feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]]
        # print 'after deletion, ', feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]]
        blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]

    return feature_all_filtered
