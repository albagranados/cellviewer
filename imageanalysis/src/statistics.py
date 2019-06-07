import math, time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# params = {
#     'text.latex.preamble': ['\\usepackage{gensymb}'],
#     'image.origin': 'lower',
#     'image.interpolation': 'nearest',
#     'image.cmap': 'gray',
#     'axes.grid': False,
#     'savefig.dpi': 150,  # to adjust notebook inline plot size
#     'axes.labelsize': 8, # fontsize for x and y labels (was 10)
#     'axes.titlesize': 8,
#     'font.size': 8, # was 10
#     'legend.fontsize': 6, # was 10
#     'xtick.labelsize': 8,
#     'ytick.labelsize': 8,
#     'text.usetex': True,
#     'figure.figsize': [3.39, 2.10],
#     'font.family': 'serif',
# }
params = {
    'axes.labelsize': 20,  # fontsize for x and y labels (was 10)
    'font.size': 24,  # was 10
    'legend.fontsize': 8, # was 10
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'font.family': 'serif',
}
mpl.rcParams.update(params)


def plot_hist(data, fig=None, ax=None, bins=None, xscale='linear', yscale={},
              xlabel='', cmap=None, num_bins=100, xticks=None,
              xtickslabel=None, alpha=1, discrete=0, xlim=None, ylim=None, color='k'):
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
    if fig is None:
        fig, ax = plt.subplots()

    if xscale is 'log':
        data = data[np.where((data != float('inf')) & (data > 0))]  # eliminate areas that are set to inf, -1 or 0
    else: data = data[np.where((data != float('inf')) & (data >= 0))]

    unique, counts = np.unique(data, return_counts=True)

    if bins is not None:
        if fig is None: fig, ax = plt.subplots(figsize=(10, 2))
        if cmap is None:
            n, bins, patches = ax.hist(data, bins=bins, histtype='bar', rwidth=1,
                                    weights=np.zeros_like(data) + 1. / data.size, color='w', ec='k')  # , color='k')
        # ax.set_aspect('equal')  #, adjustable='box')
        if cmap is not None:
            n, bins, patches = ax.hist(data, bins=bins, histtype='bar', rwidth=2,
                                    weights=np.zeros_like(data) + 1. / data.size, color='w', ec='w')  # , color='k')
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
                ax.bar(bar_xpositions, bar_heights, width=width_log, edgecolor=color, fc=(1, 1, 1, 0))
                ax.set_xscale("log")
            else:
                ax.bar(bar_xpositions, bar_heights, width=width_lin, edgecolor=color, fc=(1,1,1, 0))
        if yscale is 'log':
            ax.set_yscale('log')

    plt.ylabel(r'frequency'); ax.hold(1)
    plt.xlabel(xlabel)
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
    result (dict): cmeans, cmedians, ... see: https://mpl.org/api/_as_gen/mpl.axes.Axes.violinplot.html
    """

    if not isinstance(data, list): data = [data]  # correct if not list input data

    # alpha=0.3
    if yscale is 'log':
        data = [bp[np.where((bp != float('inf')) & (bp > 0))] for bp in data]
        values = [np.log10(bp) for bp in data]
        ylabel = ylabel + '$\quad$' + r'[$10**$]'
    else:
        values = [bp[np.where((bp != float('inf')) & (bp >= 0))] for bp in data]
    fig, ax = plt.subplots(figsize=(10, 5))  # figsize=(len(data)*1.5, len(data)*3))
    if bptype == 'violin':   # plot violin plot
        # print ylabel
        # print [np.mean(data) for data in values]
        result = ax.violinplot(values, widths=widths, showmeans=1, showmedians=True, showextrema=True)
        for ii, pc in enumerate(result['bodies']):
            pc.set_alpha(alpha)
            if cmap is None:
                pc.set_facecolor('gray'); pc.set_linewidth(1)
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
                # print 'mean in boxplot = ', result['cmeans']
                pc.set_linestyle('--')
    else:  # plot box plot
        result = ax.boxplot(values)

    # ax.set_aspect('equal', adjustable='box')

    # # adding horizontal grid lines
    # ax.yaxis.grid(True)
    # ax.set_xticks([y + 1 for y in range(len(values))])
    # ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if xticklabels == '': xticklabels=[str(i+1) for i in range(len(values))]
    # xticklabels = [r'10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%']
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
                      cluster_cluster=1, strength=1):
    """
    plot and save graphs regarding detection of clusters (diameter, densities...)

    Input:
    ----------------
    feature_all: it can be a list of features extracted from a data set or just one feature from one sample.
    """
    if not isinstance(feature_all, list): feature_all = [feature_all]  # only one sample
    if not isinstance(file_dirs, list): file_dirs = [file_dirs]  # only one sample

    # if cluster_cluster:
    #     for feature in feature_all: feature['clusterincluster'] = count_clusterincluster(feature, pixel_size)
    feature_all_dir1 = []  # will be non-empty if >1 file_dirs, i.e., if we want to plot two histograms
    # corresponding to statistics of cluster descriptors of two experiments
    for ii, file_dir in enumerate(file_dirs):
        if ii == 0: feature_all_dir0 = [feature for feature in feature_all if feature['file_dir'] == file_dir]
        elif ii == 1: feature_all_dir1 = [feature for feature in feature_all if feature['file_dir'] == file_dir]
        elif ii > 1: raise ValueError('error in statistic_descrip: only two histograms can be ovelaid.')

    blob_diameters_all_dir0 = np.asarray([diameter for feature in feature_all_dir0 for jj, diameter in enumerate(
                                feature['diameter']) if len(feature['histogram_descr'][jj]) > 0])
    cluster_density_all_dir0 = np.asarray([(feature['argmaxgrad'][0].shape[0]-feature['histogram_descr'].count([]))/
                                           float(pixel_size**2 *feature['image'].shape[0] * feature['image'].shape[1]) for feature in
                                          feature_all_dir0])
    strength_all_dir0 = np.asarray([strength for feature in feature_all_dir0
                            for jj, strength in enumerate(feature['strength']) if len(feature['histogram_descr'][jj]) > 0])
    if feature_all_dir1:
        blob_diameters_all_dir1 = np.asarray([diameter for feature in feature_all_dir1 for jj, diameter in
                                              enumerate(feature['diameter']) if len(feature['histogram_descr'][jj]) > 0])
        cluster_density_all_dir1 = np.asarray([(feature['argmaxgrad'][0].shape[0]-feature['histogram_descr'].count([]))/
                                           float(pixel_size**2 *feature['image'].shape[0] * feature['image'].shape[1]) for feature in
                                          feature_all_dir1])
        strength_all_dir1 = np.asarray([strength for feature in feature_all_dir1
                                        for jj, strength in enumerate(feature['strength']) if
                                        len(feature['histogram_descr'][jj]) > 0])
    if radius:
        fig, ax = plot_hist(0.5*blob_diameters_all_dir0, xscale='log', xlabel=r'blob radius R [nm]', discrete=1,
                            color='k', ylim=[0, 0.4])
        # plot_boxplot(0.5 * blob_diameters_all, bptype='violin', ylabel=r'blob radius R [nm]')
        if feature_all_dir1:
            plot_hist(0.5 * blob_diameters_all_dir1, fig=fig, ax=ax, xscale='log', xlabel=r'blob radius R [nm]',
                      discrete=1, color='r', ylim=[0, 0.4])
        plt.savefig(savefile_suffix + '_blobradius.pdf', bbox_inches='tight')
    if area:
        fig, ax = plot_hist(np.pi*(blob_diameters_all_dir0/2.)**2, xscale='log', xlabel=r'blob area [nm2]',
                            ylim=[0, 0.3], discrete=1)
        # plot_boxplot(np.pi * (blob_diameters_all / 2.) ** 2, bptype='violin', ylabel=r'blob area [nm2]')
        if feature_all_dir1:
            plot_hist(np.pi * (blob_diameters_all_dir1 / 2.) ** 2, fig=fig, ax=ax, xscale='log',
                      xlabel=r'blob area [nm$^2$]', discrete=1, ylim=[0, 0.3], color='r')
        plt.savefig(savefile_suffix + '_blobareas.pdf', bbox_inches='tight')
        # errorbar_featureresponse(feature, dict_sift, xlabel=r'blob diameter [nm]')
    # print '\tNumber of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0]
    # print '\tDensity of clusters:\t', feature_all[0]['argmaxgrad'][0].shape[0] / float(pixel_size *
    #        feature_all[0]['image'].shape[0] * feature_all[0]['image'].shape[1]), '[cluster/nm2]'
    if density_cluster:
        print cluster_density_all_dir0
        fig, ax = plot_hist(cluster_density_all_dir0, xscale='log', xlabel=r'$\kappa$ ['r'blobs/nm$^2$]')
        if feature_all_dir1:
            plot_hist(cluster_density_all_dir1, fig=fig, ax=ax, color='r', xscale='log',
                      xlabel=r'$\kappa$ ['r'blobs/nm$^2$]')
        plt.savefig(savefile_suffix + '_blobdensitieskappa.pdf', bbox_inches='tight')
    if strength:
        fig, ax = plot_hist(strength_all_dir0, xscale='log', xlabel=r'strength - laplacian, max$_t\{\Delta_{'
                                                                r'\gamma-norm}\}$', ylim=[0,0.1])
        if feature_all_dir1:
            plot_hist(strength_all_dir1, fig=fig, ax=ax, color='r', xscale='log',
                      xlabel=r'strength - laplacian, max$_t\{\Delta_{\gamma-norm}\}$', ylim=[0,0.1])
        plt.savefig(savefile_suffix + '_strength.pdf', bbox_inches='tight')
    if ispp:
        density_all_dir0 = np.asarray([den for feature in feature_all_dir0 for jj, den in enumerate(feature[
                                                            'density']) if len(feature['histogram_descr'][jj]) > 0])
        vor_areas_all_dir0 = np.asarray([area for feature in feature_all_dir0 for area in feature['vor'].areas])
        num_loc_all_dir0 = np.asarray([nl for feature in feature_all_dir0 for jj, nl in enumerate(feature[
                                        'number_localizations']) if len(feature['histogram_descr'][jj]) > 0])
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
            density_all_dir1 = np.asarray([den for feature in feature_all_dir1 for jj, den in enumerate(feature[
                                'density']) if len(feature['histogram_descr'][jj]) > 0])
            vor_areas_all_dir1 = np.asarray([area for feature in feature_all_dir1 for area in feature['vor'].areas])
            num_loc_all_dir1 = np.asarray([nl for feature in feature_all_dir1 for jj, nl in
                                           enumerate(feature['number_localizations']) if len(feature['histogram_descr'][jj]) > 0])
        if area_voronoi:
            fig, ax = plot_hist(vor_areas_all_dir0, xscale='log', num_bins=50, xlabel=r'Voronoi polygon area [nm$^2$]')
            # plot_boxplot(vor_areas_all, xscale='log', bptype='violin', ylabel=r'Voronoi polygon area [nm$^2$]')
            if feature_all_dir1:
                plot_hist(vor_areas_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, color='r', ylim=[0,0.3],
                          xlabel=r'Voronoi polygon 'r'area [nm$^2$]')
            plt.savefig(savefile_suffix + '_voronoiareas.pdf', bbox_inches='tight')
        if num_loc:
            fig, ax = plot_hist(num_loc_all_dir0, xscale='log', num_bins=50, xlim=[1, 10**3], ylim=[0,0.2],
                                xlabel=r'N$^{cluster}$ ['r'points/blob]')
            # plot_boxplot(num_loc_all, xscale='log', bptype='violin', ylabel=r'N$^{cluster}$ [points/cluster]')
            if feature_all_dir1:
                plot_hist(num_loc_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, xlim=[1, 10 ** 3], ylim=[0,0.2],
                          color='r', xlabel=r'N$^{cluster}$ ['r'points/blob]')
            plt.savefig(savefile_suffix + '_blobnumloc.pdf', bbox_inches='tight')
        if density:
            fig, ax = plot_hist(density_all_dir0, xscale='log', num_bins=50, ylim=[0,0.2],
                                xlabel=r'blob density [points/nm$^2$]')
            # plot_boxplot(density_all, bptype='violin', ylabel=r'cluster densities $\rho^{cluster}$ [points/nm$^2$]')
            if feature_all_dir1:
                plot_hist(density_all_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, ylim=[0,0.15], color='r',
                          xlabel=r'blob density ['r'points/nm$^2$]')
            plt.savefig(savefile_suffix + '_blobdensities.pdf', bbox_inches='tight')
        if nnd:
            dist_all_features_dir0 = np.asarray([cc for feature in feature_all_dir0 for jj, cc in enumerate(feature[
                'nnd']) if len(feature['histogram_descr'][jj]) > 0])
            # dist_all_features_dir0 = nnd_feature(feature_all_dir0, pixel_size)
            fig, ax = plot_hist(dist_all_features_dir0, xscale='log', num_bins=50, ylim=[0,0.2],
                                xlim=[0, 400], xlabel=r'nnd [nm]')
            if feature_all_dir1:
                # dist_all_features_dir1 = nnd_feature(feature_all_dir1, pixel_size)
                dist_all_features_dir1 = np.asarray([cc for feature in feature_all_dir1 for jj, cc in enumerate(feature[
                                                    'nnd']) if len(feature['histogram_descr'][jj]) > 0])
                plot_hist(dist_all_features_dir1, fig=fig, ax=ax, xscale='log', num_bins=50, color='r', ylim=[0, 0.2],
                          xlim=[0, 400], xlabel=r'nnd [nm]')
            plt.savefig(savefile_suffix + '_nnd.pdf', bbox_inches='tight')
        if cluster_cluster:
            clusterincluster_all_dir0 = np.asarray([cc for feature in feature_all_dir0 for jj, cc in enumerate(feature[
                                                   'clusterincluster']) if len(feature['histogram_descr'][jj]) > 0])
            if feature_all_dir1:
                clusterincluster_all_dir1 = np.asarray([cc for feature in feature_all_dir1 for jj, cc in enumerate(
                                            feature['clusterincluster']) if len(feature['histogram_descr'][jj]) > 0])
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
        import mpl.colors as colors
        import mpl.cm as cmx

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


def scatterplot_descriptors(feature_all, clusters, cmap=None, savefile_suffix=None,
                            filter_pos=None, cmap_charact=None, pixel_size=5, histogram=None):
    """
    This functions plots in the PC1-PC2 space the clustered distribution of all detected features to create the
    vocabulary.
    """

    from sklearn.decomposition import PCA

    if not isinstance(feature_all, list): feature_all = [feature_all]

    if histogram is None:
        histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
        histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
        if len(histogram) == 0:
            print 'error: please, compute orientation and descriptors histograms for SIFT'
            return None
    else: print '\t\tNB: input histogram should be output weighted histogram from function create_codewords.'

    # labels = clusters.labels_
    labels = np.asarray([feature['word_label'][jj][ii] for feature in feature_all for jj, orientation in
              enumerate(feature['orientation']) for ii, ori in enumerate(orientation)])
    n_clusters = np.unique(labels).shape[0]
    print '\tcomputing PCA for visualization...'
    pca = PCA(n_components=2)

    X = np.append(histogram, clusters.cluster_centers_, axis=0)  # shape (n_samples, n_features)
    X_reduced = pca.fit_transform(X)  # array-like, shape (n_samples, n_components)
    print '\t\t\tPC explained variance: ', pca.explained_variance_
    if filter_pos is not None:
        X_reduced_filtered = np.append(X_reduced[0:X_reduced.shape[0]-2][filter_pos],
                                       X_reduced[X_reduced.shape[0]-2:], axis=0)
        # X_reduced_filtered = X_reduced[filter_pos]
        histogram_filtered = histogram[filter_pos]
        labels_filtered = labels[filter_pos]
        del histogram, labels, X_reduced
        histogram = histogram_filtered; labels = labels_filtered; X_reduced = X_reduced_filtered
    fig, ax = plt.subplots(); alpha = 0.7
    if cmap_charact is None:
        if cmap is None: cmap = plt.cm.get_cmap('Set1')
        shape = ['o', 'v', 'D', '^', 's', 'P', '>', '<']  # plot marker as a function of experiemnt
        file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
        # print '\t\t\t(scatter plot: o ->', file_dirs[0].split('/')[-2], '; v ->', file_dirs[1].split('/')[-2], ')'
        labels_exp = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all
                                  for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                  enumerate(orientation)])
        for ii, fd in enumerate(file_dirs):
            pos = np.where(labels_exp == ii)[0]
            sc = plt.scatter(X_reduced[pos, 0], X_reduced[pos, 1], c=labels[pos],
                             alpha=alpha, facecolors="None", cmap=cmap, s=10, marker=shape[ii])  # cmap)
            # sc = plt.scatter(X_reduced[0:histogram.shape[0], 0], X_reduced[0:histogram.shape[0], 1], c=labels,
            #                  alpha=alpha, facecolors="None", cmap=cmap, s=10)  # cmap)  # same markers
            sc.set_facecolor('none'); plt.hold(True)
            # plt.scatter(X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 0],
            #             X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 1], marker='*', s=10 ** 2,
            #             color=[cmap(ii) for ii in range(n_cluster)])
        cmap_charact = r'1, \ldots,' + str(n_clusters)  # " ".join(str(x) for x in (np.unique(labels)+1))
    else:
        cmap = plt.cm.get_cmap('jet')
        if cmap_charact is 'strength':
            labels = np.asarray([feature[cmap_charact][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in enumerate(orientation)])
            savefile_suffix += '_strength'
        elif cmap_charact is 'area':
            labels = np.log10(np.asarray([np.pi*(0.5*feature['diameter'][ii])**2 for feature in feature_all
                                 for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                 enumerate(orientation)]))
            cmap_charact = r'blob area [nm$^2$] [10**]'; savefile_suffix += '_blobareas'
        elif cmap_charact is 'density':
            labels = np.log10(np.asarray([feature[cmap_charact][ii] for feature in feature_all
                                    for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                    enumerate(orientation)]))
            cmap_charact = r'blob density [points/nm$^2$] [10**]'; savefile_suffix += '_blobdensities'
        elif cmap_charact is 'num_loc':
            labels = np.log10(np.asarray([feature['number_localizations'][ii] for feature in feature_all
                                   for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                   enumerate(orientation)]))
            cmap_charact = r'N$^{cluster}$ ['r'points/blob] [10**]'; savefile_suffix += '_blobnumloc'
        elif cmap_charact is 'nnd':
            distances_nnd_features_all = nnd_feature(feature_all, pixel_size)
            aux = 0; nnds = []
            for feature in feature_all:
                for ii, orientation in enumerate(feature['orientation']):
                    for jj, ori in enumerate(orientation): nnds.append(distances_nnd_features_all[aux])
                    aux += 1
            labels = np.log10(np.asarray(nnds))
            cmap_charact = 'nnd [nm] [10**]'; savefile_suffix += '_nnd'
        elif cmap_charact is 'cluster_density':
            labels = np.log10(np.asarray([feature['argmaxgrad'][0].shape[0] /
                                            float(pixel_size * feature['image'].shape[0] * feature['image'].shape[1])
                                            for feature in feature_all for ii, orientation in
                                            enumerate(feature['orientation']) for jj, ori in enumerate(orientation)]))
            cmap_charact = r'$\kappa$ [blobs/nm$^2$] [10**]'; savefile_suffix += '_blobdensitieskappa'
        elif cmap_charact is 'experiment':
            file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
            # color = ['red', 'black']
            # print '\t\t\t(scatter plot:', color[0], '->', file_dirs[0].split('/')[-2], ';', color[1], '->', \
            #     file_dirs[1].split('/')[-2], ')'
            colormap = mpl.cm.viridis; normalize = mpl.colors.Normalize(vmin=0, vmax=len(file_dirs) - 1)
            color = colormap(normalize(range(len(file_dirs))))
            # for kk, col in enumerate(color):
            #     print '\t\t\tscatter plot (experiment): ', col, '->', file_dirs[kk].split('/')[-2]
            labels = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all
                                      for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                      enumerate(orientation)]); alpha = 0.5
            cmap = mpl.colors.ListedColormap(color)
            cmap_charact = ' - '.join([''.join(list(fd.split('/')[-2])[0:3]) for fd in file_dirs])
            savefile_suffix += '_experiment'
        elif cmap_charact is 'cluster_cluster':
            for feature in feature_all: feature['clusterincluster'] = count_clusterincluster(feature, pixel_size)
            labels = np.log10(np.asarray([feature['clusterincluster'][ii] for feature in feature_all
                                              for ii, orientation in enumerate(feature['orientation'])
                                              for jj, ori in enumerate(orientation)]))
            cmap_charact = r'blobs in blobs [blobs/blobs]'; savefile_suffix += '_blobsinblobs'
        sc = plt.scatter(X_reduced[0:histogram.shape[0], 0], X_reduced[0:histogram.shape[0], 1], c=labels,
                         alpha=alpha, facecolors="None", cmap=cmap, s=8)
    plt.colorbar(sc, label=cmap_charact, ticks=[])
    sc.set_facecolor('none'); plt.hold(True)
    # aux = 0
    # orientation = feature.get('orientation')
    # argmaxgrad = feature.get('argmaxgrad')
    # for i in range(argmaxgrad[0].size):
    #     for ii, ori in enumerate(np.asarray(orientation[i])):  # if empty => no in X_reduced
    #         # ax.annotate('(%d,%d), ori=%.1f' % (argmaxgrad[0][i], argmaxgrad[1][i], ori),
    #         #             xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
    #         ax.annotate('%d' % labels[aux], xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
    #         aux += 1
    # plt.plot(cluster_centers_reduced[:, 0], cluster_centers_reduced[:, 1], 'k*', markersize=10)
    # for ii in range(n_cluster):
    #     pos = np.where(siftclusters.labels_== ii)[0][0]
    #     ax.annotate('%d' % ii, xy=(X_reduced[pos, 0], X_reduced[pos, 1]), size=20)
    #     # ax.annotate('%d' % ii, xy=(X_reduced[histogram.shape[0] + ii, 0],
    #     #                            X_reduced[histogram.shape[0] + ii, 1]), size=40)
    # # ax.annotate('%.d', (X_reduced[hist_ind, s0], X_reduced[hist_ind, 1]), size=10)
    plt.xlabel(r'PC$_1$'); plt.ylabel(r'PC$_2$')
    # ax.set_ylim([-0.8, 0.8]); ax.set_xlim([-0.8, 0.6])
    ax.set_aspect('equal', adjustable='box')
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


def create_codewords(algorithm, feature_in, kwargs_sift={}, n_cluster=3, weight=0, reduced_sample=1):
    """
    This function

    Input:
    ---------
    feature_all: list of feature (argmaxgrad, orientation, scale, histogram_descr, ...) of all training images

    Output:
    ---------
    codewords_clusters, feature_codewords = dict(feature_all=feature_all_codewords, codewords_ref=codewords_ref)
    codewords_ref is an array of size equal number of centroids, of triplet containing (kk=image_no, jj=argmaxgrad, ii=orientationpos)

    """
    from sklearn.cluster import KMeans
    from sklearn.cluster import AgglomerativeClustering
    import copy

    print 'Creating vocabulary with %d words...' % n_cluster; start_time = time.time()

    if not isinstance(feature_in, list): feature_in = [feature_in]

    if reduced_sample:
        feature_all_clustering, feature_pos_reduced = reduce_dimensionality(feature_in, method='class_equalsize')  # random!
    else: feature_all_clustering = feature_in
    histogram_descr_all = [feature['histogram_descr'] for feature in feature_all_clustering]
    histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])  # some hist can be []
    if len(histogram) == 0:
        print 'error: please, compute orientation and descriptors histograms for SIFT'; return None
    if algorithm is 'kmeans':
        if weight:
            # # # experiment
            # file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
            # labels = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all
            #                      for ii, orientation in enumerate(feature['orientation']) for jj, ori in
            #                      enumerate(orientation)])
            # # blob size
            labels = np.asarray([feature['strength'][ii] for feature in feature_all_clustering
                                 for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                 enumerate(orientation)])
            # # # normalize to maximum value in histograms to avoid extra numerical error
            lmax = np.max(labels); hmax = np.max(histogram)
            labels = 1./lmax*labels*hmax
            histogram = np.append(histogram, np.array(labels).reshape(histogram.shape[0], 1),axis=1)
            # histogram /= np.linalg.norm(histogram); histogram[np.where(histogram > 0.2)] = 0.2
            # histogram /= np.linalg.norm(histogram); w1 = 1; w2 = 1
            w1 = (np.sqrt(np.mean(2*np.var(histogram[:, 0:-1], 0))))**-1
            w2 = (np.sqrt(2*np.var(histogram[:, -1], 0)))**-1
            histogram[:, 0:-1] = w1*histogram[:, 0:-1]; histogram[:, -1] = w2*histogram[:, -1]
            print '\t\tNB: modify output and return histogram_weighted.'
        codewords_clusters = KMeans(n_clusters=n_cluster, random_state=0, init="k-means++").fit(histogram)
    elif algorithm is 'hierarchical':
        codewords_clusters = AgglomerativeClustering(n_clusters=n_cluster, affinity='euclidean', linkage='ward').fit(
            histogram)

    print ("Done (time =  %.2f seconds)" % (time.time() - start_time))

    if algorithm is 'hierarchical':
        print '\t\tHierarchical clustering for codewords building. Beware that no centroids are computed.'
        codewords_clusters.cluster_centers_ = np.empty(shape=(0, 128))
        return codewords_clusters, histogram  # no data for centroids of clusters

    # # create feature list with the codewords only -------
    n_cluster = np.unique(codewords_clusters.labels_).shape[0]
    feature_all_codewords = copy.deepcopy(feature_all_clustering); codewords_ref = np.repeat(None, n_cluster)
    feature_pos_out = np.array([], dtype=int)
    for label in range(n_cluster):
        feature_pos_label = np.where(codewords_clusters.labels_ == label)
        distances = np.linalg.norm(histogram[feature_pos_label] - codewords_clusters.cluster_centers_[label], axis=1)
        if feature_pos_label[0].shape[0] > 1:
            feature_pos_label_outball = np.delete(range(feature_pos_label[0].shape[0]), np.argmin(distances))
            feature_pos_out = np.append(feature_pos_out, feature_pos_label[0][feature_pos_label_outball])
        image_no, feature_no, ori_no = feature_reference(feature_pos_label[0][np.argmin(distances)], feature_all_clustering)
        codewords_ref[label] = (image_no, feature_no, ori_no)
    feature_pos_out = np.sort(feature_pos_out)
    argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
    for kk, feature in enumerate(feature_all_clustering):
        for jj, ori in enumerate(feature['orientation']):
            for ii, o in enumerate(ori):
                image_nos.append(kk); argmaxgrad_feature_pos.append(jj); ori_feature_pos.append(ii)
    aux = 0; blob_loc_prev = None; blob_image_prev = None
    for pos in feature_pos_out:
        pos_ori = pos
        if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # feature in same blob
            aux += 1; pos_ori -= aux
        else: aux = 0
        feature_all_codewords[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
            np.delete(feature_all_codewords[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]],
                      ori_feature_pos[pos_ori])
        del feature_all_codewords[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
            pos_ori]]
        blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]
    # # -----------------------------------------------------

    return codewords_clusters, dict(feature_all=feature_all_codewords, codewords_ref=codewords_ref)


def codewords_statistics(feature_all, cmap, xlabel=r'codeword', ylabel='', pixel_size=1,
                            savefile_suffix=None, pt='_words', file_dirs=None,
                            radius=1, area=1, num_loc=1, density=1, nnd=0, strength=1, cluster_density=0, experiment=0,
                            cluster_cluster=1,
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

    labels = np.asarray([feature['word_label'][jj][ii] for feature in feature_all for jj, orientation in
                                        enumerate(feature['orientation']) for ii, ori in enumerate(orientation)])
    diameters = np.asarray([feature['diameter'][ii] for feature in feature_all
                                    for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                    enumerate(orientation)])
    feature_pos = np.where((diameters > diameter_range[0]) | (diameters < diameter_range[1]))[0]
    labels_filtered = labels[feature_pos]
    n_clusters = count_features(feature_all)
    fig, ax = plot_hist(np.asarray(labels_filtered), alpha=0.7, cmap=cmap, xlabel=xlabel,
                        # xticks=np.unique(labels_filtered),
                        bins=range(n_clusters+1),  #np.append(np.unique(labels_filtered),
                        # np.max(labels_filtered)+1)*0.4-0.2,
                        #xtickslabel=[str(e) for e in np.unique(labels_filtered)+1]
                        )
    fig.set_figheight(len(np.unique(labels_filtered))*4); fig.set_figwidth(len(np.unique(labels_filtered))*2)
    if savefile_suffix is not None: plt.savefig(savefile_suffix + '_histogram' + pt + '.pdf', bbox_inches='tight')
    print '\tStatistical test of significance between codewords characteristics:'
    if radius:
        diameters_words = [diameters[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(diameters_words) == 2:
            sign, accept = significance_test(diameters_words[0], diameters_words[1])
            print '\t\t\tRadius: H0 is', bool(accept), '-->', sign
            text = sign
        plot_boxplot([dd*0.5 for dd in diameters_words], bptype='violin', xlabel=xlabel,
                     ylabel=r'blob radius R [nm]', cmap=cmap, widths=0.5, alpha=0.7, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobradius' + pt + '.pdf', bbox_inches='tight')
    if area:
        diameters_words = [diameters[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(diameters_words) == 2:
            sign, accept = significance_test(diameters_words[0], diameters_words[1])
            print '\t\t\tArea: H0 is ', bool(accept), '-->', sign
            text = sign
        plot_boxplot([np.pi*(dd*0.5)**2 for dd in diameters_words], bptype='violin', xlabel=xlabel,
                     ylabel=r'blob area [nm$^2$]', cmap=cmap, widths=0.8, alpha=0.7, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobareas' + pt + '.pdf', bbox_inches='tight')
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
                     alpha=0.7, text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_strength' + pt + '.pdf', bbox_inches='tight')
    if density:
        densities = np.asarray([feature['density'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        densities_words = [densities[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(densities_words) == 2:
            sign, accept = significance_test(densities_words[0], densities_words[1])
            print '\t\t\tDensity: H0 is ', bool(accept), '-->', sign
            text = sign
        plot_boxplot(densities_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'blob density [points/nm$^2$]', cmap=cmap, widths=0.8, alpha=0.7, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobdensities' + pt + '.pdf', bbox_inches='tight')
    if num_loc:
        num_locs = np.asarray([feature['number_localizations'][ii] for feature in feature_all
                                for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                enumerate(orientation)])
        num_locs_words = [num_locs[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = None
        if len(num_locs_words) == 2:
            sign, accept = significance_test(num_locs_words[0], num_locs_words[1])
            text = sign
            print '\t\t\tNumber of localizations: H0 is', bool(accept), '-->', sign
        plot_boxplot(num_locs_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'N$^{cluster}$ ['r'points/blob]', cmap=cmap, widths=0.8, alpha=0.7, yscale='log',
                     text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_blobnumloc' + pt + '.pdf', bbox_inches='tight')
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
                     ylabel=r'nnd [nm]', cmap=cmap, widths=0.8, alpha=0.7, yscale='log', text=text)
        if savefile_suffix is not None: plt.savefig(savefile_suffix + '_nnd' + pt + '.pdf', bbox_inches='tight')
    if cluster_density:
        cluster_densities = np.asarray([(feature['argmaxgrad'][0].shape[0] - feature['histogram_descr'].count([]))/
                                        float(pixel_size*feature['image'].shape[0]*feature['image'].shape[1])
                                        for feature in feature_all for ii, orientation in
                                        enumerate(feature['orientation']) for jj, ori in enumerate(orientation)])
        cluster_densities_words = [cluster_densities[np.where(labels == l)[0]] for l in np.unique(labels)]
        plot_boxplot(cluster_densities_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'$\kappa$ ['r'blobs/nm$^2$]', cmap=cmap, widths=0.8, alpha=0.7, yscale='log')
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_blobdensitieskappa' + pt + '.pdf', bbox_inches='tight')
    if experiment:
        file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
        print '\t\t\t\t(NB: words statist. descrip. EXPERIMENT coded as: ', \
            [fd.split('/')[-2] for fd in file_dirs], ')'
        experiments = np.asarray([np.where(feature['file_dir'] == np.asarray(file_dirs))[0][0]
                                  for feature in feature_all
                                  for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                                  enumerate(orientation)])
        experiments_words = [experiments[np.where(labels == l)[0]] for l in np.unique(labels)]
        text = ''
        if len(experiments_words) == 2:
            sign, accept = significance_test(experiments_words[0], experiments_words[1])
            text = sign
            print '\t\t\tExperiment: H0 is', bool(accept), '-->', sign
        # plot_boxplot(experiments_words, bptype='violin', xlabel=xlabel,
        #              ylabel=file_dirs[0].split('/')[-2]+'-'+file_dirs[1].split('/')[-2], cmap=cmap, widths=0.8,
        #              alpha=0.7, text=text)
        plot_boxplot(experiments_words, bptype='violin', xlabel=xlabel,
                     ylabel=ylabel, cmap=cmap, widths=0.8, alpha=0.7, text=text)
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_experiment' + pt + '.pdf', bbox_inches='tight')
    if cluster_cluster:
        for feature in feature_all:
            feature['clusterincluster'] = count_clusterincluster(feature, pixel_size)
        cluster_in_clusters = np.asarray([feature['clusterincluster'][ii] for feature in feature_all
                                          for ii, orientation in enumerate(feature['orientation'])
                                          for jj, ori in enumerate(orientation)])
        cluster_in_clusters_words = [cluster_in_clusters[np.where(labels == l)[0]] for l in np.unique(labels)]
        if len(cluster_in_clusters_words) == 2:
            sign, accept = significance_test(cluster_in_clusters_words[0], cluster_in_clusters_words[1])
            text = sign
            print '\t\t\tBlobs in blobs: H0 is', bool(accept), '-->', sign
        plot_boxplot(cluster_in_clusters_words, bptype='violin', xlabel=xlabel,
                     ylabel=r'blobs in blobs [blobs/blobs]', yscale='log', cmap=cmap, widths=0.5, alpha=0.7, text=text)
        if savefile_suffix is not None:
            plt.savefig(savefile_suffix + '_blobsinblobs' + pt + '.pdf', bbox_inches='tight')


def sift2shape(siftclusters, n_hist=16, n_bins_descr=8, cmap=None, savefile_suffix=None, alpha=0.7):
    """
    After unsupervised-clustered algorithm for vocabulary contruction, build histograms and/or shapes of the estimated
    words (or centroids)

    Input:
    -----------------

    """

    # mpl.rcParams.update({'font.size': 10})

    if cmap is None: cmap = plt.cm.get_cmap('Set1')
    # mpl.rcParams.update({'font.size': 10})

    for no, hist in enumerate(siftclusters.cluster_centers_):
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


def filter_features(criterion, feature_all, kwargs, siftclusters=None):
    """
    Function to filter detected features considering PC1-PC2 plots with different criteria

    criterion == 'histogram_distance' refers to distance to cluster centers in the PC space. This requires parameter
    ball_radius = percentile all distnaces in that label
    criterion == 'blob_size' filters according to blob diameter. This requires parameter diameter,
    criterion == 'experiment' filters according to experiment (file_dir). This require parameter experiment_name
    criterion == 'strength'
    criterion == 'numlocs' with minimum number of localizations per blob as input parameter
    criterion == 'words', with word_id as input parameter
    """
    import copy
    if not isinstance(feature_all, list): feature_all = [feature_all]

    diameters = np.asarray([feature['diameter'][ii] for feature in feature_all
                 for ii, orientation in enumerate(feature['orientation']) for jj, ori in enumerate(orientation)])
    num_features = diameters.shape[0]
    print '\t\t\tnumber of features before filtering: ', num_features
    siftclusters_filtered = None
    if criterion is 'histogram_distance':
        ball_radius_percentile = kwargs.get('ball_radius_percentile', 50)
        if siftclusters is None: print 'error: need siftclusters computations.'
        labels = siftclusters.labels_
        cluster_centers = siftclusters.cluster_centers_  # siftclusters.cluster_centers_[0] = centroid of cluster 0
        histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
        histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
        if len(histogram) == 0:
            print 'error: please, compute orientation and descriptors histograms for SIFT'; return None
        feature_pos_out = np.array([], dtype=int)
        for label in np.unique(labels):
            feature_pos_label = np.where(labels == label)
            distances = np.linalg.norm(histogram[feature_pos_label] - cluster_centers[label], axis=1)
            ball_radius = np.percentile(distances, ball_radius_percentile)
            print '\t\t\tFiltering: cluster label ', label, ', filter PC at a distance to centeroid less than ', \
                ball_radius
            feature_pos_label_outball = np.where(distances > ball_radius)
            feature_pos_out = np.append(feature_pos_out, feature_pos_label[0][feature_pos_label_outball[0]])
        feature_pos_out = np.sort(feature_pos_out)
        siftclusters_filtered = copy.deepcopy(siftclusters)
        siftclusters_filtered.labels_ = np.delete(siftclusters_filtered.labels_, feature_pos_out)
    if criterion is 'words':
        word_id = kwargs.get('word_id', 0)
        if siftclusters is None: print 'error: need siftclusters computations.'
        labels = siftclusters.labels_
        feature_pos_out = np.where(labels != word_id)[0]
        siftclusters_filtered = copy.deepcopy(siftclusters)
        siftclusters_filtered.labels_ = np.delete(siftclusters_filtered.labels_, feature_pos_out)
    elif criterion is 'blob_size':
        diameter_range = kwargs.get('diameter_range', [-float('inf'), float('inf')])
        feature_pos_out = np.where((diameters < diameter_range[0]) | (diameters > diameter_range[1]))[0]
    elif criterion is 'experiment':
        experiment_name = kwargs.get('experiment_name', 0)
        feature_pos_out = np.where(np.asarray([feature['file_dir'] != experiment_name for kk, feature in enumerate(feature_all)
                                   for orientation in feature['orientation'] for ori in orientation], dtype=int) == 1)[0]
    elif criterion is 'strength':
        threshold_percent = kwargs.get('threshold_percent', None)
        threshold_value = kwargs.get('threshold_value', None)
        # strengths = np.asarray([feature['strength'][ii] for feature in feature_all
        #                         for ii, orientation in enumerate(feature['orientation']) for jj, ori in
        #                         enumerate(orientation)])
        # if threshold_percent is not None:
        #     threshold = np.max(strengths) - threshold_percent * (np.max(strengths) - np.min(strengths))
        # elif threshold_value is not None: threshold = threshold_value
        # feature_pos_out = np.where(strengths < threshold)[0]
        strengths_imno = np.asarray([(feature['strength'][ii], imno) for imno, feature in enumerate(feature_all)
                                        for ii, orientation in enumerate(feature['orientation'])
                                        for jj, ori in enumerate(orientation)])
        for imno, feature in enumerate(feature_all):
            strengths = strengths_imno[np.where(strengths_imno[:, 1] == imno)[0]][:, 0]
            # fig, ax = plot_hist(strengths, xscale='log', xlabel=r'strength - laplacian, max$_t\{\Delta_{'
            #                                                             r'\gamma-norm}\}$', ylim=[0, 0.1])
            if threshold_percent is not None:
                threshold = np.max(strengths) - threshold_percent * (np.max(strengths) - np.min(strengths))
            elif threshold_value is not None: threshold = threshold_value
            else:  # if no threshold is given, filter by clustering strength distribution - remove lower strength
                from sklearn.cluster import KMeans
                from sklearn.decomposition import PCA
                clusters = KMeans(n_clusters=2, random_state=0, init="k-means++"). \
                    fit(np.reshape(strengths, (strengths.shape[0], 1)))
                labels = clusters.labels_
                position_out = np.where(np.argmin(clusters.cluster_centers_) == labels)[0]
                threshold = np.max(strengths[position_out])
            if imno == 0:
                feature_pos_out = np.where((strengths_imno[:, 0] < threshold) & (strengths_imno[:, 1] == imno))[0]
            else:
                feature_pos_out = np.concatenate([feature_pos_out, np.where((strengths_imno[:, 0] < threshold) &
                                                                            (strengths_imno[:, 1] == imno))[0]])
            print '\t\t\t\t filtering feature in image ', imno, 'with strength threshold ', threshold
            # # filter by clustering strength distribution image to image
            # from sklearn.cluster import KMeans
            # from sklearn.decomposition import PCA
            # clusters = KMeans(n_clusters=2, random_state=0, init="k-means++").\
            #                     fit(np.reshape(strengths,(strengths.shape[0],1)))
            # labels = clusters.labels_
            # position_out = np.where(np.argmin(clusters.cluster_centers_) == labels)[0]
            # print '\t\t\t\t filtering feature in image ', imno, 'with strength threshold ', np.max(strengths[
            #                                                                                             position_out])
            # if imno == 0: feature_pos_out = np.where(strengths_imno[:, 1] == imno)[0][position_out]
            # else:
            #     feature_pos_out = np.concatenate([feature_pos_out,
            #                                       np.where(strengths_imno[:, 1] == imno)[0][position_out]])
    elif criterion is 'numloc':
        min_numloc = kwargs.get('min_numloc', 1)
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
        if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # feature in same blob
            aux += 1; pos_ori -= aux
        else: aux = 0
        feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
            np.delete(feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]],
                      ori_feature_pos[pos_ori])
        del feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
            pos_ori]]
        blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]

    print '\t\t\tnumber of features after filtering: ', num_features-feature_pos_out.shape[0]

    return feature_all_filtered, feature_pos_filtered, siftclusters_filtered


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


def count_clusterincluster(feature, pixel_size, max_clusterincluster=50):
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


def compute_silhouette(data, labels, cmap=None, savefile_suffix=None):

    from sklearn.metrics import silhouette_samples, silhouette_score

    print '\tComputing average silhouette_score... '
    silhouette_avg = silhouette_score(data, labels)
    sample_silhouette_values = silhouette_samples(data, labels)
    # silhouette_avg = np.median(sample_silhouette_values)
    n_clusters = np.unique(labels).shape[0]
    print '\tDone. The average silhouette_score is :', silhouette_avg
    if savefile_suffix is not None:
        # fig, ax = plt.subplots();
        fig, ax = plt.subplots(figsize=(n_clusters * 1.5, n_clusters * 3))
        ax.set_xlim([0, len(data) + (n_clusters + 1) * 10])
        x_lower = 10; xticks_pos = []
        cluster_max = []  # if we are interested in percentil representing
                          # average among maxima of each cluster (some/all/none above average silhouette)
        for i in range(n_clusters):
            ith_cluster_silhouette_values = sample_silhouette_values[labels == i]
            ith_cluster_silhouette_values.sort()

            cluster_max.append(max(ith_cluster_silhouette_values))

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            x_upper = x_lower + size_cluster_i

            ax.fill_between(np.arange(x_lower, x_upper), 0, ith_cluster_silhouette_values,
                              facecolor=cmap(i), edgecolor=cmap(i), alpha=0.7)
            # ax.text(x_lower + 0.5 * size_cluster_i, -0.05, str(i+1))
            xticks_pos.append(x_lower + 0.5 * size_cluster_i)
            x_lower = x_upper + 10
        cluster_max.sort()
        silhouette_percentile = round((len([i for i, x in enumerate(cluster_max) if x < silhouette_avg]) / float(len(
            cluster_max))), 2)

        print '\t\t(', silhouette_percentile*100, '% of clusters fall below this value)'

        # ax.set_ylabel(); ax.set_xlabel()

        # The vertical line for average silhouette score of all the values
        ax.axhline(y=silhouette_avg, color='k', linestyle="--")
        plt.setp(ax, xticks=xticks_pos, xticklabels=[str(i + 1) for i in range(n_clusters)],
                 ylabel=r'silhouette coefficient $s(i)$', xlabel="cluster - visual word")

        # ax.set_xticks([])  # Clear the yaxis labels / ticks
        fig.savefig(savefile_suffix + '_silhouette.pdf', bbox_inches='tight')

    return silhouette_avg, silhouette_percentile


def reduce_dimensionality(feature_all, n_cluster=50, method='kmeans', pca_comp=30):

    ini_time = time.time()
    import copy
    from sklearn.cluster import MiniBatchKMeans
    if not isinstance(feature_all, list): feature_all = [feature_all]

    print '\t\treducing dimensionality of training set by', method, '...'
    if method is ('minibatch' or 'kmeans'):
        print method
        if method is 'minibatch':
            histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
            histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
            siftclusters = MiniBatchKMeans(init='k-means++', n_clusters=n_cluster, batch_size=50,
                                           n_init=10, max_no_improvement=10, verbose=0).fit(histogram)
        elif method is 'kmeans':
            siftclusters, histogram = create_codewords('kmeans', feature_all, n_cluster=n_cluster)
            print '\t\t\t done.'
        labels = siftclusters.labels_
        cluster_centers = siftclusters.cluster_centers_
        feature_pos_out = np.array([], dtype=int)
        num_features = count_features(feature_all)
        print '\t\t\tnumber of features before filtering: ', num_features
        for label in np.unique(labels):
            print '\n\n\n\n label = ', label
            feature_pos_label = np.where(labels == label)
            distances = np.linalg.norm(histogram[feature_pos_label] - cluster_centers[label], axis=1)  # if only one
            if feature_pos_label[0].shape[0] > 1:
                feature_pos_label_outball = np.delete(range(feature_pos_label[0].shape[0]), np.argmin(distances))
                feature_pos_out = np.append(feature_pos_out, feature_pos_label[0][feature_pos_label_outball])
        feature_pos_out = np.sort(feature_pos_out)
        feature_pos_filtered = np.delete(range(num_features), feature_pos_out)
        feature_all_filtered = copy.deepcopy(feature_all)
        argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
        for ii, feature in enumerate(feature_all):  # for all images/windows
            for jj, ori in enumerate(feature['orientation']):
                for kk, o in enumerate(ori):
                    image_nos.append(ii)
                    argmaxgrad_feature_pos.append(jj)
                    ori_feature_pos.append(kk)
        aux = 0; blob_loc_prev = None; blob_image_prev = None
        for pos in feature_pos_out:
            pos_ori = pos
            if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # feature in same blob
                aux += 1; pos_ori -= aux
            else: aux = 0
            feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
                np.delete(feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]],
                          ori_feature_pos[pos_ori])
            del feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
                pos_ori]]
            blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]

        print '\t\t\tnumber of features after filtering: ', num_features-feature_pos_out.shape[0]
        return feature_all_filtered, feature_pos_filtered

    elif method is 'pca':  # see (Fernando et al., 2012), p. 901
        from sklearn.decomposition import PCA
        histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
        histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])
        pca = PCA(n_components=pca_comp)
        histogram_reduced = pca.fit_transform(histogram)  # array-like, shape (n_samples, n_components)
        print '\t\t\tPC explained variance with', pca_comp, 'components is: ', np.sum(pca.explained_variance_)
        feature_all_filtered = copy.deepcopy(feature_all)
        num_feat = 0
        for feat_no, feature in enumerate(feature_all):
            for keypoint_pos, histogram in enumerate(feature['histogram_descr']):
                for keypoint_descr, hist_sub in enumerate(histogram):
                    feature_all_filtered[feat_no]['histogram_descr'][keypoint_pos][keypoint_descr]=histogram_reduced[
                        num_feat]
                    num_feat += 1
        feature_pos_filtered = range(num_feat)

        return feature_all_filtered, feature_pos_filtered

    elif method is 'class_equalsize':
        import random
        file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
        labels = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all
                             for ii, orientation in enumerate(feature['orientation']) for jj, ori in
                             enumerate(orientation)])
        class_size = []  # vector with number of features per class/label-as-above
        for l in np.unique(labels): class_size.append(np.where(labels == l)[0].shape[0])
        equalsize = min(class_size)
        print '\t\t\tSelect', equalsize, 'features per class.'
        num_features = sum(class_size)
        feature_pos_out = np.array([], dtype=int)
        for l in np.unique(labels):
            print '\t\t\texperiment/class', l, '(', file_dirs[l].split('/')[-2], ') is reduced by', \
                (class_size[l]-equalsize), 'features.'
            feature_pos_label_out = random.sample(np.where(labels == l)[0], k=(class_size[l]-equalsize))  # random!
            if len(feature_pos_label_out) > 0:
                feature_pos_out = np.append(feature_pos_out, feature_pos_label_out)
        feature_pos_out = np.sort(feature_pos_out)
        feature_pos_filtered = np.delete(range(num_features), feature_pos_out)
        feature_all_filtered = copy.deepcopy(feature_all)
        argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
        for ii, feature in enumerate(feature_all):  # for all images/windows
            for jj, ori in enumerate(feature['orientation']):
                for kk, o in enumerate(ori):
                    image_nos.append(ii); argmaxgrad_feature_pos.append(jj); ori_feature_pos.append(kk)
        aux = 0; blob_loc_prev = None; blob_image_prev = None
        for pos in feature_pos_out:
            pos_ori = pos
            if blob_loc_prev == argmaxgrad_feature_pos[pos] and blob_image_prev == image_nos[pos]:  # feature in same blob
                aux += 1; pos_ori -= aux
            else: aux = 0
            feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]] = \
                np.delete(feature_all_filtered[image_nos[pos]]['orientation'][argmaxgrad_feature_pos[pos]],
                          ori_feature_pos[pos_ori])
            del feature_all_filtered[image_nos[pos]]['histogram_descr'][argmaxgrad_feature_pos[pos]][ori_feature_pos[
                pos_ori]]
            blob_loc_prev = argmaxgrad_feature_pos[pos]; blob_image_prev = image_nos[pos]

        return feature_all_filtered, feature_pos_filtered

    print ("\t\tDone (total time =  %.2f seconds)" % (time.time() - ini_time))


def quantize_descriptors(feature_all, clusters, assignment='hard'):

    """
    Image representation.

    This function generates an array-histogram for each image-roi (in feature['codewords_hist'] ), used as input for the
    classifier. Requires a previously computed vocabulary.

    Also new attribut feature['word_label'].

    Creates new entries in the list of dictionaries feature_all as feature['codewords_hist]. Input is a mutable
    object - a list of dictionaries

    hard assignment (Fernando et al., 2012)
    """

    import copy

    histogram_descr_all = [feature['histogram_descr'] for feature in feature_all]
    histogram = np.asarray([hist_sub for histogram_descr in histogram_descr_all for hist in histogram_descr for hist_sub in hist])  # some hist can be []
    labels = clusters.predict(histogram)   # array with labels corresponding to each feature vector in each image
    n_clusters = np.unique(labels).shape[0]

    print 'Quantize image features into a', n_clusters, '-dimensional vector...'; ini_time = time.time()
    if assignment is 'hard':
        f0 = 0; n_features = 0  # f0-to-fn is image (roi) range
        for kk, feature in enumerate(feature_all):
            feature['word_label'] = []
            #  for jj, orientations in enumerate(feature['orientation']):
                #  feature['word_label'][jj] = []
                #  for ii, ori in enumerate(orientations): feature['word_label'][jj].append(labels[n_features+ii])
            for jj, orientations in enumerate(feature['orientation']):
                feature['word_label'].append(labels[(f0+n_features):(f0+n_features+len(orientations))])
                n_features += len(orientations)
            fn = f0 + n_features; words = labels[f0:fn].tolist()
            feature['codewords_hist'] = np.asarray([words.count(l) for l in range(n_clusters)])
            f0 = fn; n_features = 0
    print ("Done (time =  %.2f seconds)" % (time.time() - ini_time))


def count_features(feature_all):

    """
    This function counts the number of features (with no-empty orientation)
    """

    if not isinstance(feature_all, list): feature_all = [feature_all]

    num_features = np.asarray([kk for kk, feature in enumerate(feature_all) for ii, orientation in
                                enumerate(feature['orientation']) for jj, ori in enumerate(orientation)]).shape[0]
    return num_features


def plot_codewords(feature_all, feature_codewords, plot_num_features=3,
                   savefile_suffix=None, dict_inputfile={}, dict_image={}, dict_sift={}):

    """
    This function plots image patches corresponding to the closes keypoint to the pre-computed cluster centroids in create_codewords (output feature_codewords)

    The codewords are sorted in descending order according to the size of its membership (Fei2005).
    The appearance of 10 codewords selected from the top 20 most likely codewords for this category model.

    """
    from src import iprocessing as iproc
    from src import vprocessing as vproc
    from src import utilities as util
    from scipy.spatial import Voronoi

    print 'Plot codewords histograms and most likely codewords for each class/category/experiment... '; ini_time = \
        time.time()
    labels = np.asarray([feature['word_label'][jj][ii] for feature in feature_all for jj, orientation in
                         enumerate(feature['orientation']) for ii, ori in enumerate(orientation)])
    n_clusters = np.unique(labels).shape[0]
    file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
    # experiments = np.asarray([np.where(feature['file_dir'] == np.asarray(file_dirs))[0][0]
    #                           for feature in feature_all
    #                           for ii, orientation in enumerate(feature['orientation']) for jj, ori in
    #                           enumerate(orientation)])
    codewords_ref = feature_codewords['codewords_ref']; feature_all_codewords = feature_codewords['feature_all']
    codewords_hist_exp = np.zeros(shape=(file_dirs.shape[0], n_clusters)); n_features_class = np.zeros(file_dirs.shape[0])
    for cl, fd in enumerate(file_dirs):
        n_features = 0
        for feature in [feat for feat in feature_all if feat['file_dir'] == fd]:
            codewords_hist_exp[cl, :] += feature['codewords_hist']
            n_features += count_features(feature)
        codewords_hist_exp[cl, :] = 1. / n_features * codewords_hist_exp[cl, :]
        n_features_class[cl] = n_features
    if savefile_suffix is not None:
        ylim = [0, 1.1 * np.max(codewords_hist_exp)]
        for cl, fd in enumerate(file_dirs):
            colormap = mpl.cm.viridis; normalize = mpl.colors.Normalize(vmin=0, vmax=len(file_dirs) - 1)
            color = colormap(normalize(range(len(file_dirs))))
            cmap = mpl.colors.ListedColormap(np.asarray([color[np.where(feature_all[image_no]['file_dir'] ==
                             file_dirs)[0][0]] for (image_no, feature_no, ori_no) in codewords_ref]))
            # cmap_charact = [f.split('/')[-2].split('_')[0] for f in file_dirs] # print cmap_charact
            fig, ax = plot_hist(np.asarray([el for l, rep in
                                            enumerate(n_features_class[cl] * codewords_hist_exp[cl, :]) for el
                                            in np.repeat(l, rep)]), bins=range(n_clusters + 1),
                                cmap=cmap, xlabel=r'codewords',  # xticks=np.unique(labels)*0.2,alpha=0.7,
                                xtickslabel=None,  # [str(e) for e in np.unique(labels) + 1],
                                xlim=[-1, n_clusters], ylim=ylim, alpha=1)
            plt.savefig(savefile_suffix + '_expercm' + '_' + fd.split('/')[-2] + '_codewords_hist.pdf', bbox_inches='tight')

            cmap = util.discrete_cmap(n_clusters, "jet")
            fig, ax = plot_hist(np.asarray([el for l, rep in
                                            enumerate(n_features_class[cl] * codewords_hist_exp[cl, :]) for el
                                            in np.repeat(l, rep)]), bins=range(n_clusters+1),
                                            cmap=cmap, xlabel=r'codewords',  # xticks=np.unique(labels)*0.2,alpha=0.7,
                                            xtickslabel=None,  # [str(e) for e in np.unique(labels) + 1],
                                            xlim=[-1, n_clusters], ylim=ylim, alpha=1)
            plt.savefig(savefile_suffix + '_' + fd.split('/')[-2] + '_codewords_hist.pdf', bbox_inches='tight')
            #
            # centroids_ids = np.argsort(codewords_hist_exp[cl, :])[-plot_num_features:][::-1]  # order codewords from
            # # less frequent to most. Take indexes of last largest.
            # print '\n\tClass =', fd.split('/')[-2], '\n\tFrequent codewords labels =', centroids_ids+1
            # plt.figure(figsize=(8, 4))
            # for cw, (image_no, feature_no, ori_no) in enumerate(codewords_ref[centroids_ids]):
            #     print '\t\tPlot frequent codeword (image_no, feature_no, ori_no)=', (image_no, feature_no, ori_no)
            #     feature_image = feature_all[image_no]
            #     # image=feature_image['image']; vor=feature_image['vor']
            #     dict_inputfile['file_dir'] = feature_image['file_dir']
            #     dict_inputfile['file_name'] = feature_image['file_name']
            #     if dict_inputfile.get('ispp'):
            #         data = util.pointpattern(); data.read(dict_inputfile, plot=0)
            #         points = data.points; vor = Voronoi(points)
            #     vproc.compute_parameters(vor, dict_inputfile)  # new attribute in vor object: vor.areas
            #     image = vproc.densities_interpolate(vor, scale_pixel_size=dict_image.get('scale_pixel_size'),
            #                                         interpolate_method=dict_image.get('interpolate_method'),
            #                                         fill_value=0.0, density_transform=dict_image['detect_densitytransform'])
            #     bx = feature_image['argmaxgrad'][0][feature_no]; by = feature_image['argmaxgrad'][1][feature_no]
            #     sc = feature_image['scale'][feature_no]
            #     plt.subplot(1, centroids_ids.shape[0], cw + 1)
            #     plt.imshow(image[int(bx-2.5*np.sqrt(sc)):int(bx+2.5*np.sqrt(sc)), int(by-2.5*np.sqrt(sc)):int(by+2.5*np.sqrt(sc))].T,
            #                cmap='jet', interpolation='none'); plt.xticks(()); plt.yticks(())
            # plt.subplots_adjust(0.08, 0.02, 0.92, 0.85, 0.08, 0.23); plt.show()
            # plt.savefig(savefile_suffix + '_' + fd.split('/')[-2] + '_codewords.pdf', bbox_inches='tight')

    print ("Done (time =  %.2f seconds)" % (time.time() - ini_time))

    return 1  # codewords_hist_exp


def feature_reference(index, feature_all):

    """
    For a given position-index in a labelfeature array (such as clusters.labels_), this function returns the
    corresponding (image/roi number (kk), argmax position (jj), orientation position (ii)
    """

    argmaxgrad_feature_pos = []; image_nos = []; ori_feature_pos = []
    for kk, feature in enumerate(feature_all):
        for jj, ori in enumerate(feature['orientation']):
            for ii, o in enumerate(ori):
                image_nos.append(kk); argmaxgrad_feature_pos.append(jj); ori_feature_pos.append(ii)

    return image_nos[index], argmaxgrad_feature_pos[index], ori_feature_pos[index]


def classify(feature_all, feature_all_test):

    """"
    array X of size [n_samples, n_features] holding the training samples, and an array y of class labels
    (strings or integers), size [n_samples]:
    """

    from sklearn import svm

    print 'Test images classification...'; ini_time = time.time()
    file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
    n_classes = file_dirs.shape[0]
    n_samples = len(feature_all); n_features = feature_all[0]['codewords_hist'].shape[0]
    y = []; X = []
    print '\tClasses codes: ', [(cl, fd.split('/')[-2]) for cl, fd in enumerate(file_dirs)]
    print 'Known classes of test images:'
    for cl, fd in enumerate(file_dirs):
        for feature in [feat for feat in feature_all if feat['file_dir'] == fd]:
            X.append(1./count_features(feature)*feature['codewords_hist'])
            y.append(cl)
            print feature['file_dir'].split('/')[-2], y[-1]
    clf = svm.SVC(gamma='auto')
    clf.fit(np.asarray(X).reshape(n_samples, n_features), np.asarray(y))

    test = []; n_samples_test = len(feature_all_test)
    for cl, fd in enumerate(file_dirs):
        for feature in [feat for feat in feature_all_test if feat['file_dir'] == fd]:
            test.append(1./count_features(feature)*feature['codewords_hist'])
    labels_predict = clf.predict(np.asarray(test).reshape(n_samples_test, n_features))

    print 'Prediction of test images:'
    for no, feature in enumerate(feature_all_test):
        print feature['file_dir'].split('/')[-2], labels_predict[no]

    print ("Done (time =  %.2f seconds)" % (time.time() - ini_time))

    return 1


def scatterplot_bowvector(feature_all, feature_all_test=None, savefile_suffix=None):
    """
    This functions plots in the PC1-PC2 space the clustered distribution of all detected features to create the
    vocabulary.
    """

    from sklearn.decomposition import PCA
    from src import utilities as util

    if not isinstance(feature_all, list): feature_all = [feature_all]

    histograms = np.asarray([1./count_features(feature)*feature['codewords_hist'] for feature in feature_all])

    pca = PCA(n_components=2)
    X_reduced = pca.fit_transform(histograms)  # array-like, shape (n_samples, n_components)

    file_dirs = np.unique(np.asarray([feature['file_dir'] for feature in feature_all]))
    print '\tClasses codes: ', [(cl, fd.split('/')[-2]) for cl, fd in enumerate(file_dirs)]
    labels = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all])

    fig, ax = plt.subplots()
    colormap = mpl.cm.viridis; normalize = mpl.colors.Normalize(vmin=0, vmax=len(file_dirs) - 1)
    color = colormap(normalize(range(len(file_dirs)))); cmap = mpl.colors.ListedColormap(color)
    sc = plt.scatter(X_reduced[:, 0], X_reduced[:, 1], c=labels, alpha=0.5, cmap=cmap, s=8)
    sc.set_facecolor('none')

    plt.hold(True)
    histograms_test = np.asarray([1./count_features(feature)*feature['codewords_hist'] for feature in feature_all_test])
    pca = PCA(n_components=2); X_reduced_test = pca.fit_transform(histograms_test)
    labels_test = np.asarray([np.where(feature['file_dir'] == file_dirs)[0][0] for feature in feature_all_test])
    sc = plt.scatter(X_reduced_test[:, 0], X_reduced_test[:, 1], c=labels_test, alpha=0.5, cmap=cmap, s=10,
                     marker='^')
    sc.set_edgecolor('none')

    plt.xlabel(r'PC$_1$'); plt.ylabel(r'PC$_2$'); plt.xticks()
    # ax.set_ylim([-0.8, 0.8]); ax.set_xlim([-0.8, 0.6])
    ax.axes.get_xaxis().set_ticks([]); ax.axes.get_yaxis().set_ticks([])
    ax.set_aspect('equal', adjustable='box')

    cbar = util.colorbar(sc)
    cbar.ax.set_yticklabels('')
    cbar.set_label(r' - '.join([''.join(list(fd.split('/')[-2])[0:3]) for fd in file_dirs]))

    if savefile_suffix is not None:
        plt.savefig(savefile_suffix + '_scplot_bowvector2' + '.pdf', bbox_inches='tight')
