import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# configure latex plots
matplotlib.rcParams['text.usetex'] = True  # uselatex labels
matplotlib.rcParams['font.family'] = 'serif'


def plot_hist(data, hist_scale='log', num_bins=100, bins=None):
    """
    This function plots histograms. Range: min-max of data.

    Input:
    ----------------
    hist_scale (string): log or lin xaxis of the histogram. Default = 'log'
            i.e., bins=np.logspace(np.log10(ini), np.log10(np.max(data)), num=num_bins)
                  bins=np.linspace(ini, np.max(data), num=num_bins)
    num_bins (float): number of bins. Default = 100
    bins (1d array): if it is not None, then bins=bins and 'lin' scale

    """
    import numpy as np
    import matplotlib.pyplot as plt

    ini = np.min(data)

    if bins is not None:
        plt.figure(), plt.hist(data, bins=bins, histtype='step',
                               weights=np.zeros_like(data) + 1. / data.size, color='k')
        # plt.gca().set_xscale("log")
        plt.ylabel(r'frequency')
        return

    if hist_scale is 'log':
        if ini <= 0: ini = 1
        plt.figure(), plt.hist(data, bins=np.logspace(np.log10(ini), np.log10(np.max(data)), num=num_bins),
                               histtype='step', weights=np.zeros_like(data) + 1. / data.size, color='k')
        # plt.xlim(1, 10**3)
        plt.gca().set_xscale("log")
        plt.ylabel(r'frequency')
    else:
        if ini < 0: ini = 0
        plt.figure(), plt.hist(data, bins=np.linspace(ini, np.max(data), num=num_bins), histtype='step',
                               weights=np.zeros_like(data) + 1. / data.size, color='k')
        plt.ylabel(r'frequency')


def sample_statistics(data):

    mean = round(np.average(data), 2)
    std = round(np.std(data, ddof=1), 2)  # sum()/(N-ddof)
    median = round(np.median(data), 2)

    return {'mean': mean, 'std': std, 'median': median}


def errorbar_featureresponse(feature, dict_analysis):
    """
    Plot the error bar strength response vs blob diameter

    Input:
    -------------------
    feature (dictionary)
    dict_analysis (dictionary)

    """

    import math

    featurestrength = feature.get('featurestrength')
    argmaxgrad = feature.get('argmaxgrad')
    tnew = feature.get('tnew')
    scale_range = feature.get('scale_range')

    err = np.zeros(shape=scale_range.shape)
    y = np.zeros(shape=scale_range.shape)

    for ii in range(len(scale_range)):
        pos_equal_scale = np.where(tnew == scale_range[ii])
        strength = featurestrength[argmaxgrad[0][pos_equal_scale], argmaxgrad[1][pos_equal_scale]]
        y[ii] = np.median(strength)
        if math.isnan(y[ii]): y[ii] = 0
        err[ii] = np.std(strength)

    plt.figure(); plt.errorbar(dict_analysis.get('analysis_pixel_size')*4*np.sqrt(scale_range), y, err, fmt='ko')
    plt.xlabel(r'blob diameter [nm]'); plt.ylabel('max$_t\{\Delta_{\gamma-norm}\}$')

    return y


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
    photons = pc*area  #  number of photons for each molecule in the molecule list. The photon conversion factor,
                        # eg. 0.41 for STORM and 0.14 for NSTORM
    N = photons # total photon count, number of photons in the image
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
            ax = plt.plot(points[ii][0], points[ii][1], color=deviation, markersize=1)
            plt.hold(True)
        plt.colorbar(scalarMap, label='$\sqrt{2\mu_x}$')
        plt.hold(False)

    return rms


def siftdescr_analysis(feature, kwargs_sift={}, n_cluster=3, init="k-means++", max_iter=300, n_jobs=-1,
                       compute_pca=True, plot_graphics=True, fig={}, ax={}):
    """
    This function

    Input:
    ---------
    feature dictionary with sift elements, e.g., 'histogram_descr'.

    Output:
    ---------
    """

    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA

    orientation = feature.get('orientation')
    argmaxgrad = feature.get('argmaxgrad')
    histogram_descr = feature.get('histogram_descr')
    histogram = np.asarray([hist for allhist_feature in histogram_descr for hist in allhist_feature])
    if len(histogram) == 0:
        print 'warning: please, compute orientation and descriptors histograms for SIFT'
        return None

    n_bins_descr = kwargs_sift.get('n_bins_descr', 4)
    n_hist = kwargs_sift.get('n_hist', 4)

    kmeans = KMeans(n_clusters=n_cluster, random_state=0, init=init).fit(histogram)
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_

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

    if compute_pca:
        print 'computing PCA for visualization (n_clusters = %.d)...' % n_cluster
        pca = PCA(n_components=2)
        X = np.append(histogram, cluster_centers, axis=0)  #  shape (n_samples, n_features)
        X_reduced = pca.fit_transform(X)  # array-like, shape (n_samples, n_components)
        # histogram_reduced = pca.fit_transform(histogram)  # array-like, shape (n_samples, n_components)
        # cluster_centers_reduced = pca.fit_transform(cluster_centers)  # array-like, shape (n_samples, n_components)

        if plot_graphics:
            if fig is None: fig, ax = plt.subplots()
            ax.scatter(X_reduced[0:histogram.shape[0], 0], X_reduced[0:histogram.shape[0], 1], c=labels,
                       alpha=0.8, cmap=plt.cm.get_cmap('Set1'))
            fig.hold(True)
            aux = 0
            for i in range(argmaxgrad[0].size):
                for ii, ori in enumerate(np.asarray(orientation[i])):
                    # ax.annotate('(%d,%d), ori=%.1f' % (argmaxgrad[0][i], argmaxgrad[1][i], ori),
                    #             xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
                    ax.annotate('(%d,%d)' % (argmaxgrad[0][i], argmaxgrad[1][i]),
                                xy=(X_reduced[aux, 0], X_reduced[aux, 1]))
                    aux += 1
            ax.plot(X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 0],
                    X_reduced[histogram.shape[0]:histogram.shape[0] + n_cluster, 1], 'k*', markersize=10)
            # # ax.plot(cluster_centers_reduced[:,0], cluster_centers_reduced[:, 1],
            # #            'k*', markersize=10)
            plt.xlabel(r'PCA$_1$'); plt.ylabel(r'PCA$_2$')
            plt.show(); fig.hold(False)

            # arg = 0; hist_ind = 0
            # for hist in histogram_descr:
            #     if not hist:  # skip that blob, out of range
            #         arg += 1
            #         continue
            #     bx = argmaxgrad[0][arg]; by = argmaxgrad[1][arg]
            #     t = tnew[arg]
            #     arg += 1
            #     for jj in range(len(hist)):
            #         ax.annotate('%.0f,%.0f' % (bx, by), (X_reduced[hist_ind, 0]-0.01, X_reduced[hist_ind, 1]-0.01), size=10)
            #         ax.annotate('%.0f' % t, (X_reduced[hist_ind, 0], X_reduced[hist_ind, 1]), size=10)
            #         hist_ind += 1

    return kmeans

#  kmeans.predict([[0, 0], [4, 4]])

