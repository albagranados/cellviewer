import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams["text.usetex"] = True; matplotlib.rcParams['font.family'] = 'serif'  # configure latex plots
matplotlib.rcParams.update({'font.size': 16})


def plot_hist(data, fig=None, ax=None, bins=None, xscale='linear', yscale={},
              xlabel={}, cmap=None, num_bins=100, xticks=None,
              xtickslabel=None, alpha=1, discrete=0, xlim=None, ylim=[0, 0.95], color='k'):
    """
    This function plots histograms. Range: min-max of data. Compared to plot_hist, this new version computes the
    histogram with bar function and considering the x-range discrete/categorical if necessary (scale is typically
    discrete=1, e.g., blob radius. But 'continuous'/no discrete in e.g., number of points/blob).
    The width in x-log scale is fixed.

    Input:
    ----------------
    discrete (boolean): if true or 1, then the x-range of histogram is computed with bar function to avoid plotting
                        non-existing data
    scale (string): log or lin xaxis of the histogram. Default = 'log'
            i.e., bins=np.logspace(np.log10(ini), np.log10(np.max(data)), num=num_bins)
                  bins=np.linspace(ini, np.max(data), num=num_bins)
    num_bins (float): number of bins. Default = 100
    bins (1d array): if it is not None, then bins=bins and 'lin' scale

    """
    if fig is None: fig, ax = plt.subplots()

    if xscale is 'log':
        data = data[np.where((data != float('inf')) & (data > 0))]  # eliminate areas that are set to inf, -1 or 0
    else: data = data[np.where((data != float('inf')) & (data >= 0))]

    unique, counts = np.unique(data, return_counts=True)

    if bins is not None:
        n, bins, patches = ax.hist(data, bins=bins, histtype='bar', rwidth=0.4,
                                    weights=np.zeros_like(data) + 1. / data.size, color='w')  # , color='k')
        ax.set_aspect('equal')  #, adjustable='box')
        if cmap is not None:
            for c, p in zip(range(len(bins)-1), patches):
                plt.setp(p, 'facecolor', cmap(c), alpha=alpha)
        if xticks is not None:
            if xtickslabel is None: xtickslabel = xticks
            plt.xticks(xticks, xtickslabel)
        fig.show()
    else:
        if not discrete:
            ini = np.min(data)
            if xscale is 'log':
                ax.hist(data, bins=np.logspace(np.log10(ini * 0.99), np.log10(np.max(data * 1.01)), num=num_bins),
                        histtype='step', weights=np.zeros_like(data) + 1. / data.size,
                        color=color, alpha=alpha)
                ax.set_xscale("log")
            else:
                ax.hist(data, bins=np.linspace(ini, np.max(data), num=num_bins), histtype='step',
                        weights=np.zeros_like(data) + 1. / data.size, color=color, alpha=alpha)
        else:
            ini = np.min(data)
            bar_heights = 1. / np.sum(counts) * counts
            bar_xpositions = unique
            width_lin = 0.02
            if xscale is 'log':
                width_log = 10 ** width_lin * bar_xpositions - bar_xpositions
                # fig, ax = plt.subplots()
                ax.bar(bar_xpositions, bar_heights, width=width_log, edgecolor=color, fc=(1,1,1, 0))
                ax.set_xscale("log")
            else:
                ax.bar(bar_xpositions, bar_heights, width=width_lin, edgecolor=color, fc=(1,1,1, 0))
        if yscale is 'log':
            ax.set_yscale('log')

    plt.ylabel(r'frequency'); plt.xlabel(xlabel); ax.hold(1)
    ax.set_ylim(ylim); ax.set_xlim(xlim)

    return fig, ax


def plot_hist_v0(data, bins=None, xscale='lin', xlabel={}, cmap=None, num_bins=100, xticks=None, xtickslabel=None,
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

    if xscale is 'log':
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
        if xscale is 'log':
            ax.hist(data, bins=np.logspace(np.log10(ini*0.99), np.log10(np.max(data*1.01)), num=num_bins),
                                   histtype='step', weights=np.zeros_like(data) + 1. / data.size,
                                   color='k', alpha=alpha)
            ax.set_xscale("log")
        else:
            ax.hist(data, bins=np.linspace(ini, np.max(data), num=num_bins), histtype='step',
                                   weights=np.zeros_like(data) + 1. / data.size, color='k', alpha=alpha)
            # ax.hist(data, bins=num_bins, histtype='step',
            #                        weights=np.zeros_like(data) + 1. / data.size, color='k')
    plt.ylabel(r'frequency'); plt.xlabel(xlabel); ax.hold(0)
    ax.set_ylim(0, 0.9)
    # ax.hist(densities, bins='rice', histtype='step',  color='k'); plt.ylabel(r'counts')

    return fig, ax


def plot_boxplot(data, yscale='lin', bptype='violin', xticklabels='', xlabel='', ylabel='values', cmap=None,
                 widths=0.95, alpha=0.5, text=None):
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

    if yscale is 'log':
        data = [bp[np.where((bp != float('inf')) & (bp > 0))] for bp in data]
        values = [np.log10(bp) for bp in data]
        ylabel = ylabel + '$\quad$' + r'[$10**$]'
    else:
        values = [bp[np.where((bp != float('inf')) & (bp >= 0))] for bp in data]
    fig, ax = plt.subplots(figsize=(len(data)*1.5, len(data)*3))
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
                # pc.set_edgecolor(cmap(ii))
                # pc.set_facecolor('None')
                # pc.set_linewidth(2)
        for pc in [result[ii] for ii in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans')]:
            pc.set_edgecolor('black')
            pc.set_linewidth(2)
            pc.set_alpha(alpha)
            if pc == result['cmeans']:
                pc.set_linestyle('--')
    else:  # plot box plot
        result = ax.boxplot(values)

    # ax.set_aspect('equal', adjustable='box')

    # # adding horizontal grid lines
    # ax.yaxis.grid(True)
    # ax.set_xticks([y + 1 for y in range(len(values))])
    # ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if xticklabels == '': xticklabels=[str(i+1) for i in range(len(values))]
    plt.setp(ax, xticks=[y + 1 for y in range(len(values))], xticklabels=xticklabels, ylabel=ylabel, xlabel=xlabel)
    if text is not None: plt.title(text)
    plt.show()

    return result


def sample_statistics(data):

    mean = round(np.nanmean(data), 2)
    std = round(np.nanstd(data, ddof=1), 2)  # sum()/(N-ddof)
    median = round(np.nanmedian(data), 2)

    return {'mean': mean, 'std': std, 'median': median}


def statistic_descrip(feature_all, file_dirs, ispp=1, pixel_size=1, savefile_suffix=None,
                      radius=1, area=1, density_cluster=1, area_voronoi=1, num_loc=1, density=1, nnd=0,
                      cluster_cluster=1):
    """
    plot and save graphs regarding detection of clusters (diameter, densities...)

    Input:
    ----------------
    feature_all: it can be a list of features extracted from a data set or just one feature from one sample.
    """
    if not isinstance(feature_all, list): feature_all = [feature_all]  # only one sample
    if not isinstance(file_dirs, list): file_dirs = [file_dirs]  # only one sample

    if cluster_cluster:
        for feature in feature_all: feature['clusterincluster'] = count_clusterincluster(feature, pixel_size)
    feature_all_dir1 = []  # will be non-empty if >1 file_dirs, i.e., if we want to plot two histograms
    # corresponding to statistics of cluster descriptors of two experiments
    for ii, file_dir in enumerate(file_dirs):
        if ii == 0: feature_all_dir0 = [feature for feature in feature_all if feature['file_dir'] == file_dir]
        elif ii == 1: feature_all_dir1 = [feature for feature in feature_all if feature['file_dir'] == file_dir]
        elif ii > 1: raise ValueError('error in statistic_descrip: only two histograms can be ovelaid.')

    blob_diameters_all_dir0 = np.asarray([diameter for feature in feature_all_dir0 for diameter in feature['diameter']])
    cluster_density_all_dir0 = np.asarray([feature['argmaxgrad'][0].shape[0] / float(pixel_size**2 *
                                          feature['image'].shape[0] * feature['image'].shape[1]) for feature in
                                          feature_all_dir0])
    if feature_all_dir1:
        blob_diameters_all_dir1 = np.asarray([diameter for feature in feature_all_dir1 for diameter in feature['diameter']])
        cluster_density_all_dir1 = np.asarray([feature['argmaxgrad'][0].shape[0] /
                                               float(pixel_size ** 2 *feature['image'].shape[0] *feature['image'].shape[1]) for
                                               feature in feature_all_dir1])
    if radius:
        fig, ax = plot_hist(0.5*blob_diameters_all_dir0, xscale='log', xlabel=r'blob radius R [nm]', discrete=1,
                            color='k', yscale='log')
        # plot_boxplot(0.5 * blob_diameters_all, bptype='violin', ylabel=r'blob radius R [nm]')
        if feature_all_dir1:
            plot_hist(0.5 * blob_diameters_all_dir1, fig=fig, ax=ax, xscale='log', xlabel=r'blob radius R [nm]',
                      discrete=1, color='r', yscale='log')
        plt.savefig(savefile_suffix + '_blobradius.pdf', bbox_inches='tight')
    if area:
        fig, ax = plot_hist(np.pi*(blob_diameters_all_dir0/2.)**2, xscale='log', xlabel=r'blob area [nm2]',
                  yscale='log', discrete=1)
        # plot_boxplot(np.pi * (blob_diameters_all / 2.) ** 2, bptype='violin', ylabel=r'blob area [nm2]')
        if feature_all_dir1:
            plot_hist(np.pi * (blob_diameters_all_dir1 / 2.) ** 2, fig=fig, ax=ax, xscale='log',
                      xlabel=r'blob area [nm$^2$]', discrete=1, color='r', yscale='log')
        plt.savefig(savefile_suffix + '_blobareas.pdf', bbox_inches='tight')
        # errorbar_featureresponse(feature, dict_sift, xlabel=r'blob diameter [nm]')
    # print '\tNumber of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0]
    # print '\tDensity of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0] / float(pixel_size *
    #        feature_all[0]['image'].shape[0] * feature_all[0]['image'].shape[1]), '[cluster/nm2]'
    if density_cluster:
        fig, ax = plot_hist(cluster_density_all_dir0, xscale='log', xlabel=r'$\kappa$ ['r'blobs/nm$^2$]')
        if feature_all_dir1:
            plot_hist(cluster_density_all_dir1, fig=fig, ax=ax, color='r', xscale='log',
                      xlabel=r'$\kappa$ ['r'blobs/nm$^2$]')
        plt.savefig(savefile_suffix + '_blobdensitieskappa.pdf', bbox_inches='tight')
    if ispp:
        density_all_dir0 = np.asarray([den for feature in feature_all_dir0 for den in feature['density']])
        vor_areas_all_dir0 = np.asarray([area for feature in feature_all_dir0 for area in feature['vor'].areas])
        num_loc_all_dir0 = np.asarray([nl for feature in feature_all_dir0 for nl in feature['number_localizations']])
        # percentage_number_localizations = np.asarray([numloc / points.shape[0] * 100 for feature in feature_all for
                                            # numloc in feature['number_localizations']])
        # vareas_statistics = sample_statistics(vor.areas[vor.areas < float('inf')]*
        #                                       dict_inputfile['original_pixel_size']**2)
        # densities_statistics = sample_statistics(density_all)
        # feature_statistics = sample_statistics(feature['number_localizations'])
        # print '\tPercentage of localizations in clusters:\t', np.sum(percentage_number_localizations), '%'
        # print '\tVoronoi polygon areas: \t', vareas_statistics, '[nm2]'
        # print '\tNumber of loc. per blob: \t', feature_statistics, '[loc/blob]'
        # print '\tLocal blob densities: \t', densities_statistics, ' [loc/nm2]'
        if feature_all_dir1:
            density_all_dir1 = np.asarray([den for feature in feature_all_dir1 for den in feature['density']])
            vor_areas_all_dir1 = np.asarray([area for feature in feature_all_dir1 for area in feature['vor'].areas])
            num_loc_all_dir1 = np.asarray(
                [nl for feature in feature_all_dir1 for nl in feature['number_localizations']])
        if area_voronoi:
            fig, ax = plot_hist(vor_areas_all_dir0, xscale='log', num_bins=50, xlabel=r'Voronoi polygon area [nm$^2$]')
            # plot_boxplot(vor_areas_all, xscale='log', bptype='violin', ylabel=r'Voronoi polygon area [nm$^2$]')
            if feature_all_dir1:
                plot_hist(vor_areas_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, color='r', ylim=[0,0.3],
                          xlabel=r'Voronoi polygon 'r'area [nm$^2$]')
            plt.savefig(savefile_suffix + '_voronoiareas.pdf', bbox_inches='tight')
        if num_loc:
            fig, ax = plot_hist(num_loc_all_dir0, xscale='log', num_bins=50, xlim=[1, 10**3], ylim=[0,0.3],
                      xlabel=r'N$^{cluster}$ ['r'points/blob]', yscale='log')
            # plot_boxplot(num_loc_all, xscale='log', bptype='violin', ylabel=r'N$^{cluster}$ [points/cluster]')
            if feature_all_dir1:
                plot_hist(num_loc_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, xlim=[1, 10 ** 3], ylim=[0, 0.3],
                          color='r', xlabel=r'N$^{cluster}$ ['r'points/blob]', yscale='log')
            plt.savefig(savefile_suffix + '_blobnumloc.pdf', bbox_inches='tight')
        if density:
            fig, ax = plot_hist(density_all_dir0, xscale='log', num_bins=50, ylim=[0,0.3],
                                xlabel=r'blob density [points/nm$^2$]')
            # plot_boxplot(density_all, bptype='violin', ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
            if feature_all_dir1:
                plot_hist(density_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, ylim=[0,0.3], color='r',
                          xlabel=r'blob density ['r'points/nm$^2$]')
            plt.savefig(savefile_suffix + '_blobdensities.pdf', bbox_inches='tight')
        if nnd:
            dist_all_features_dir0 = nnd_feature(feature_all_dir0, pixel_size)
            fig, ax = plot_hist(dist_all_features_dir0, xscale='log', num_bins=50, ylim=[0,0.3],
                                xlim=[0, 400], xlabel=r'nnd [nm]')
            if feature_all_dir1:
                dist_all_features_dir1 = nnd_feature(feature_all_dir1, pixel_size)
                plot_hist(dist_all_features_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, color='r', ylim=[0, 0.3],
                          xlim=[0, 400], xlabel=r'nnd [nm]')
            plt.savefig(savefile_suffix + '_nnd.pdf', bbox_inches='tight')
        if cluster_cluster:
            clusterincluster_all_dir0 = np.asarray([cc for feature in feature_all_dir0 for cc in feature[
                                                   'clusterincluster']])
            if feature_all_dir1:
                clusterincluster_all_dir1 = np.asarray([cc for feature in feature_all_dir1 for cc in feature[
                                                       'clusterincluster']])
            fig, ax = plot_hist(clusterincluster_all_dir0, xscale='log', num_bins=50, xlim=[1,10**2], ylim=[0, 0.4],
                                xlabel=r'blobs in blobs [blobs/blobs]')
            if feature_all_dir1:
                plot_hist(clusterincluster_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, xlim=[1,10**2],
                          ylim=[0, 0.4], color='r', xlabel=r'blobs in blobs [blobs/blobs]')
            plt.savefig(savefile_suffix + '_blobsinblobs.pdf', bbox_inches='tight')


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


def scatterplot_vocabulary(feature_all, kmeans, n_cluster=3, cmap=None, savefile_suffix=None, filter_pos=None):
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
    if filter_pos is not None:
        X_reduced_filtered = np.append(X_reduced[0:X_reduced.shape[0]-2][filter_pos],
                                       X_reduced[X_reduced.shape[0]-2:], axis=0)
        # X_reduced_filtered = X_reduced[filter_pos]
        histogram_filtered = histogram[filter_pos]
        labels_filtered = labels[filter_pos]
        del histogram, labels, X_reduced
        histogram = histogram_filtered; labels = labels_filtered; X_reduced = X_reduced_filtered

    fig, ax = plt.subplots()
    # plt.figure()
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
    ax.set_ylim([-0.8,0.8]); ax.set_xlim([-0.8, 0.8])
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


def words_statistic_descrip(feature_all, labels, cmap, xlabel=r'cluster - visual word', pixel_size=1,
                            savefile_suffix=None,
                            radius=1, area=1, num_loc=1, density=1, nnd=0, strength=1, cluster_density=0, experiment=0,
                            cluster_cluster=0,
                            diameter_range=[-float('inf'), float('inf')]):
    """
    This function plots histograms and violin plots of the scales, densities, etc., corresponding to each cluster in the
    vocabulary

    Input:
    ---------------------

    Output:
    ---------------------

    """
    if not isinstance(feature_all, list): feature_all = [feature_all]

    diameters = np.asarray([feature['diameter'][ii] for feature in feature_all
                                    for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                    enumerate(orientation)])
    feature_pos = np.where((diameters > diameter_range[0]) | (diameters < diameter_range[1]))[0]
    labels_filtered = labels[feature_pos]
    fig, ax = plot_hist(np.asarray(labels_filtered)*0.4,
              bins=np.append(np.unique(labels_filtered), np.max(labels_filtered)+1)*0.4-0.2,
              cmap=cmap, xlabel=xlabel, xticks=np.unique(labels_filtered)*0.4,
              xtickslabel=[str(e) for e in np.unique(labels_filtered)+1], alpha=0.5)
    fig.set_figheight(len(np.unique(labels_filtered))*4); fig.set_figwidth(len(np.unique(labels_filtered))*2)
    if savefile_suffix is not None: plt.savefig(savefile_suffix + '_histogram_words.pdf', bbox_inches='tight')
    print '\tStatistical test of significance between cluster-words characteristics:'
    if radius:
        diameters_words = [diameters[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(diameters_words)==2:
            sign, accept = significance_test(diameters_words[0], diameters_words[1])
            print '\t\t\tRadius: H0 is', bool(accept), '-->', sign
            text = sign
        plot_boxplot([dd*0.5 for dd in diameters_words], bptype='violin', xlabel=xlabel,
                     ylabel=r'blob radius R [nm]', cmap=cmap, widths=0.5, alpha=0.5, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobradius_words.pdf', bbox_inches='tight')
    if area:
        diameters_words = [diameters[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(diameters_words)==2:
            sign, accept = significance_test(diameters_words[0], diameters_words[1])
            print '\t\t\tArea: H0 is ', bool(accept), '-->', sign
            text = sign
        plot_boxplot([np.pi*(dd*0.5)**2 for dd in diameters_words], bptype='violin', xlabel=xlabel,
                     ylabel=r'blob area [nm$^2$]', cmap=cmap, widths=0.8, alpha=0.5, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobareas_words.pdf', bbox_inches='tight')
    if strength:
        strengths = np.asarray([feature['strength'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        strengths_words = [strengths[np.where(labels == l)[0]] for l in np.unique(labels)]
        if len(strengths_words)==2:
            sign, accept = significance_test(strengths_words[0], strengths_words[1])
            print '\t\t\tStrength: H0 is', bool(accept), '-->', sign
            text = sign
        plot_boxplot(strengths_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'strength - laplacian, max$_t\{\Delta_{\gamma-norm}\}$', cmap=cmap, widths=0.5,
                     alpha=0.5, text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_strength_words.pdf', bbox_inches='tight')
    if density:
        densities = np.asarray([feature['density'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        densities_words = [densities[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(densities_words)==2:
            sign, accept = significance_test(densities_words[0], densities_words[1])
            print '\t\t\tDensity: H0 is ', bool(accept), '-->', sign
            text = sign
        plot_boxplot(densities_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'blob density [points/nm$^2$]', cmap=cmap, widths=0.8, alpha=0.5, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobdensities_words.pdf', bbox_inches='tight')
    if num_loc:
        num_locs = np.asarray([feature['number_localizations'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        num_locs_words = [num_locs[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(num_locs_words)==2:
            sign, accept = significance_test(num_locs_words[0], num_locs_words[1])
            text = sign
            print '\t\t\tNumber of localizations: H0 is', bool(accept), '-->', sign
        plot_boxplot(num_locs_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'N$^{cluster}$ ['r'points/blob]', cmap=cmap, widths=0.8, alpha=0.5, yscale='log',
                     text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobnumloc_words.pdf', bbox_inches='tight')
    if nnd:
        distances_nnd_features_all = nnd_feature(feature_all, pixel_size)
        aux = 0; nnds = []
        for feature in feature_all:
            for ii, orientation in enumerate(feature['orientation']):
                for jj, ori in enumerate(orientation):
                    nnds.append(distances_nnd_features_all[aux])
                aux += 1
        nnds = np.asarray(nnds)
        nnds_words = [nnds[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(nnds_words)==2:
            sign, accept = significance_test(nnds_words[0], nnds_words[1])
            text = sign
            print '\t\t\tNearest-neighbor distance: H0 is', bool(accept), '-->', sign
        plot_boxplot(nnds_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'nnd [nm]', cmap=cmap, widths=0.8, alpha=0.5, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_nnd_words.pdf', bbox_inches='tight')
    if cluster_density:
        cluster_densities = np.asarray([feature['argmaxgrad'][0].shape[0]/
                                        float(pixel_size*feature['image'].shape[0]*feature['image'].shape[1])
                                        for feature in feature_all for ii, orientation in
                                        enumerate(feature['orientation']) for jj, ori in enumerate(orientation)])
        cluster_densities_words = [cluster_densities[np.where(labels == l)[0]] for l in np.unique(labels)]
        plot_boxplot(cluster_densities_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'$\kappa$ ['r'blobs/nm$^2$]', cmap=cmap, widths=0.8, alpha=0.5, yscale='log')
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_blobdensitieskappa_words.pdf', bbox_inches='tight')
    if experiment:
        file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
        experiments = np.asarray([np.where(feature['file_dir']==file_dirs)[0][0] for feature in feature_all
                                  for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                  enumerate(orientation)])
        experiments_words = [experiments[np.where(labels == l)[0]] for l in np.unique(labels)]
        if len(experiments_words)==2:
            sign, accept = significance_test(experiments_words[0], experiments_words[1])
            text = sign
            print '\t\t\tExperiment: H0 is', bool(accept), '-->', sign
        plot_boxplot(experiments_words, bptype='violin', xlabel=xlabel,
                     ylabel=file_dirs[0].split('/')[-2]+'-'+file_dirs[1].split('/')[-2], cmap=cmap, widths=0.8,
                     alpha=0.5, text=text)
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_experiment_words.pdf', bbox_inches='tight')
    if cluster_cluster:
        for feature in feature_all:
            feature['clusterincluster'] = count_clusterincluster(feature, pixel_size)
        cluster_in_clusters = np.asarray([feature['clusterincluster'][ii] for feature in feature_all
                                          for ii, orientation in enumerate(feature['orientation'])
                                          for jj, ori in enumerate(orientation)])
        cluster_in_clusters_words = [cluster_in_clusters[np.where(labels == l)[0]] for l in np.unique(labels)]
        if len(cluster_in_clusters_words)==2:
            sign, accept = significance_test(cluster_in_clusters_words[0], cluster_in_clusters_words[1])
            text = sign
            print '\t\t\tBlobs in blobs: H0 is', bool(accept), '-->', sign
        plot_boxplot(cluster_in_clusters_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'blobs in blobs [blobs/blobs]', yscale='log', cmap=cmap, widths=0.5, alpha=0.5, text=text)
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_blobsinblobs_words.pdf', bbox_inches='tight')


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
        fig, axes = plt.subplots(nrows=4, ncols=4, sharex='col', sharey='row') # fig.suptitle(r'word %d' % (no+1))
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


def filter_features(criterion, feature_all, kmeans=None, ball_radius_percentile=50,
                    diameter_range=[-float('inf'), float('inf')], experiment_name=None, threshold_percent=0.5,
                    word_id=0, min_numloc=1):
    """
    Function to filter detected features considering PC1-PC2 plots with different criteria

    criterion == 'histogram_distance' refers to distance to cluster centers in the PC space. This requires parameter
    ball_radius = percentile all distnaces in that label
    criterion == 'blob_size' filters according to blob diameter. This requires parameter diameter,
    criterion == 'experiment' filters according to experiment (file_dir). This require parameter experiment_name
    criterion == 'strength
    criterion == 'numlocs' with minimum number of localizations per blob as input parameter
    """
    import copy

    if not isinstance(feature_all, list): feature_all = [feature_all]

    diameters = np.asarray([feature['diameter'][ii] for feature in feature_all
                 for ii, orientation in enumerate(feature['orientation']) for jj, ori in enumerate(orientation)])
    num_features = diameters.shape[0]
    kmeans_filtered = None
    if criterion is 'histogram_distance':
        if kmeans is None:
            print 'error: need kmeans computations.'
        labels = kmeans.labels_
        cluster_centers = kmeans.cluster_centers_  # kmeans.cluster_centers_[0] = centroid of cluster 0
        histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
        histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
        if len(histogram) == 0:
            print 'error: please, compute orientation and descriptors histograms for SIFT'; return None
        feature_pos_out = np.array([], dtype=int)
        for label in np.unique(labels):
            feature_pos_label = np.where(labels == label)
            distances = np.linalg.norm(histogram[feature_pos_label] - cluster_centers[label], axis=1)
            ball_radius = np.percentile(distances, ball_radius_percentile)
            print 'For cluster label ', label, ', plot ploint PC1-PC2 at a distance to cluster centers less than ', \
                ball_radius
            feature_pos_label_outball = np.where(distances > ball_radius)
            feature_pos_out = np.append(feature_pos_out, feature_pos_label[0][feature_pos_label_outball[0]])
        feature_pos_out = np.sort(feature_pos_out)
        kmeans_filtered = copy.deepcopy(kmeans)
        kmeans_filtered.labels_ = np.delete(kmeans_filtered.labels_, feature_pos_out)
    if criterion is 'words':
        if kmeans is None:
            print 'error: need kmeans computations.'
        labels = kmeans.labels_
        feature_pos_out = np.where(labels != word_id)[0]
        kmeans_filtered = copy.deepcopy(kmeans)
        kmeans_filtered.labels_ = np.delete(kmeans_filtered.labels_, feature_pos_out)
    elif criterion is 'blob_size':
        feature_pos_out = np.where((diameters < diameter_range[0]) | (diameters > diameter_range[1]))[0]
    elif criterion is 'experiment':
        feature_pos_out = np.where(np.asarray([feature['file_dir'] != experiment_name for kk, feature in enumerate(feature_all)
                                   for orientation in feature['orientation'] for ori in orientation], dtype=int) == 1)[0]
    elif criterion is 'strength':
        strengths = np.asarray([feature['strength'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        threshold = np.max(strengths) - threshold_percent * (np.max(strengths) - np.min(strengths))
        feature_pos_out = np.where(strengths < threshold)[0]
    elif criterion is 'numloc':
        num_locs = np.asarray([feature['number_localizations'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        feature_pos_out = np.where(num_locs < min_numloc)[0]
    feature_pos_filtered = np.delete(range(num_features), feature_pos_out)
    feature_all_filtered = copy.deepcopy(feature_all)
    argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
    for ii, feature in enumerate(feature_all):  # for all images/windows
        for jj, ori in enumerate(feature['orientation']):
            for kk, o in enumerate(ori):  # only count if !=[]
                image_nos.append(ii)  # for each label/feature/ori position, what is the corresponding image
                argmaxgrad_feature_pos.append(jj) # for each feature position, position of argmaxgrad in that image
                ori_feature_pos.append(kk)  # for each feature position, position of ori in that orientation
    aux = 0; blob_loc_prev = None; blob_image_prev = None
    for pos in feature_pos_out:
        pos_ori = pos
        # print 'labels[pos] = ', labels[pos]
        # print 'argmaxgrad_feature_pos[pos] = ', argmaxgrad_feature_pos[pos]
        # print 'feature_all[image_nos[pos]][orientation][argmaxgrad_feature_pos[pos]]', feature_all[image_nos[pos]][
        #     'orientation'][argmaxgrad_feature_pos[pos]]
        # print 'ori_feature_pos[pos_ori]', ori_feature_pos[pos_ori]
        # print 'ori_feature_pos[pos_ori]=', ori_feature_pos[pos_ori], '; argmaxgrad_feature_pos[pos]=', argmaxgrad_feature_pos[pos]
        if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # feature in same blob
            aux += 1; pos_ori -= aux
            # print 'IN. ori_feature_pos[pos_ori] = ', ori_feature_pos[pos_ori]
        else: aux = 0
        # print 'before deletion, feature_all_filtered: ', feature_all_filtered[image_nos[pos]]['orientation'][
        #     argmaxgrad_feature_pos[pos]]
        feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
            np.delete(feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]],
                      ori_feature_pos[pos_ori])
        del feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
            pos_ori]]
        # print 'image_nos[pos]=', image_nos[pos], 'argmaxgrad_feature_pos[pos]=', argmaxgrad_feature_pos[pos]
        # print 'after deletion, ', feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]]
        # print 'after deletion, ', feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]]
        blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]

    return feature_all_filtered, feature_pos_filtered, kmeans_filtered


def nnd_feature(feature_all, pixel_size):

    from sklearn.neighbors import NearestNeighbors

    if not isinstance(feature_all, list): feature_all = [feature_all]

    argmaxgrad_x0 = np.concatenate(np.asarray([feature['argmaxgrad'][0] for feature in feature_all]), axis=0)
    argmaxgrad_y0 = np.concatenate(np.asarray([feature['argmaxgrad'][1] for feature in feature_all]), axis=0)

    blobs_xy = np.array([argmaxgrad_x0, argmaxgrad_y0]).T
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(blobs_xy)
    distances, indices = nbrs.kneighbors(blobs_xy)

    dist_all_features = distances[:, 1] * pixel_size

    return dist_all_features


def significance_test(data1, data2, alpha=0.05, alpha1=0.01, alpha2=0.001):
    """
    The purpose of a one-way analysis of variance (one-way ANOVA) is to compare the means of two or more groups
    (the independent variable) on one dependent variable to see if the group means are significantly different from
    each other.
    """

    from scipy.stats import f_oneway

    stat, p = f_oneway(data1, data2)
    # print('Statistics=%.3f, p=%.3f' % (stat, p))
    if np.isnan(p):
        print('Statistics=%.3f, p=%.3f' % (stat, p)), ' --> exit test'
        accept = 1
        return '', accept
    if p > alpha:
        # print('Same distributions (fail to reject H0)')
        accept = 1
        return '', accept
    else:
        accept = 0
        if p < alpha2: return '***', accept
        elif p < alpha1: return '**', accept
        else: return '*', accept


def count_clusterincluster(feature, pixel_size, max_clusterincluster=100):
    """
    Count number of features within intensity-dependent and feature-based clusters (default: blobs)

    Input:
    -------------

    Output:
    -------------
    number_localizations (1d-array): following the order of argmaxgrad, each position contains the number of
    localizations within that feature (e.g., blob with radius 2*\sqrt{t})
    """

    from sklearn.neighbors import NearestNeighbors

    radius = 0.5*feature['diameter']
    blobs_xy = np.array([np.asarray(feature['argmaxgrad'][0]), np.asarray(feature['argmaxgrad'][1])]).T
    nbrs = NearestNeighbors(n_neighbors=max_clusterincluster, algorithm='ball_tree').fit(blobs_xy)
    distances, indices = nbrs.kneighbors(blobs_xy)
    distances_nm = distances * pixel_size

    # for ii in range(indices.shape[0]):
    #     np.where((radius[indices[ii, 0]] - distances_nm[ii, 1:]) > radius[indices[ii, 1:]])
    clusterincluster = np.zeros(shape=blobs_xy.shape[0])
    for ii in range(indices.shape[0]):
        cluster_in_pos = np.where((distances_nm[ii, :] < radius[indices[ii, 0]]) &
                       (radius[indices[ii, :]] < radius[indices[ii, 0]]))
        # print ii, 'at', blobs_xy[ii]
        # print radius[indices[ii,0]]
        # print distances_nm[ii, :]
        # print radius[indices[ii,:]]
        # print cluster_in_pos
        # print cluster_in_pos[0].shape[0]
        clusterincluster[ii] = cluster_in_pos[0].shape[0]

    return clusterincluster
