import os, matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# params = {
#     'axes.labelsize': 20,  # fontsize for x and y labels (was 10)
#     'font.size': 24,  # was 10
#     'legend.fontsize': 10, # was 10
#     'xtick.labelsize': 10,
#     'ytick.labelsize': 10,
#     'text.usetex': True,
#     'font.family': 'serif',
# }
# mpl.rcParams.update(params)
mpl.rcParams["text.usetex"] = True; mpl.rcParams['font.family'] = 'serif'  # configure latex plots



class pointpattern():

    #     """
    #     This function reads point pattern(s) stored in file_dir
    #
    #     Input:
    #     ----------------------------------
    #     file_dir (string): file(s) directory
    #     file_name (string): if not None, then read multple data sets in file_dir
    #     fileExt (string): extension of the pointpattern files
    #     storm (boolean): if the datafile is in storm format (precision, coordinates...)
    #     save_output_dir (string): if not None, then save the point pattern(s) in a txt files as a nx2 matrix in the output
    #                                 directory (string)
    #
    #     Output:
    #     ----------------------------------
    #     points: points (if multiple data sets, the last one) stored in a nx2 matrix
    #     points1, points2 (if more than 1 channel): stored in a nx2 matrix

    def read(self, dict_inputfile, save_output_dir=None, plot=1):

        file_dir = dict_inputfile['file_dir']; file_name = dict_inputfile['file_name']
        fileExt = dict_inputfile['fileExt']
        is_storm = dict_inputfile['is_storm']; out_channel = dict_inputfile['out_channel']

        print 'Reading dataset...'
        if file_name is not None:
            file_names = [file_name]
        else:
            file_names = os.listdir(file_dir)

        for ii in range(len(file_names)):
            file_name = file_names[ii].split(fileExt)[0]  # '2017-12-14_HeLa_DMSO_000_list_m1000DC_CTS_0851-1731_allChs'
            print '\tfile directory: "', file_dir, '"'
            print '\tfile name: "', file_name, '"'

            path = os.path.join(file_dir, file_name + fileExt)  # if we want to run in the Python Console

            if is_storm:
                self.points = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=(3, 4))  # export x_c and y_c
                print '\tAll channels in .points.'
                self.channels = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=(0,))  # export x_c and y_c
                # self.area = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=6)  # export x_c and y_c
                # self.width = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=7)  # export x_c and y_c
                # self.aspect_ratio = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=9)
                # self.bg = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=10)
                self.frame_no = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=12)

                if len(np.unique(self.channels)) > 2:
                    print '\tSTORM file with more than 2 channels. Please, double check.'
                if out_channel is not 'all':
                    for jj, ch in enumerate(out_channel):
                        if jj == 0:
                           self.points1 = self.points[np.where(self.channels == ch)]
                           print '\tChannel %d in .points1.' % ch
                           if save_output_dir is not None:
                               np.savetxt(save_output_dir + file_name + '_Ch' + str(ch) + fileExt,
                                          self.points[np.where(self.channels == out_channel[0])])
                               print '\tFile saved in *_Ch' + str(ch)
                        if jj == 1:
                            print 'in out_channel==2'
                            self.points2 = self.points[np.where(self.channels == ch)]
                            print '\tChannel %d in .points2.' % ch
                            if save_output_dir is not None:
                                np.savetxt(save_output_dir + file_name + '_Ch' + str(ch) + fileExt,
                                          self.points[np.where(self.channels == out_channel[1])])
                                print '\tFile saved in *_Ch' + str(ch)
                else:
                   if save_output_dir is True:
                       np.savetxt(save_output_dir + file_name + fileExt, self.points)
            if not is_storm:
                self.points = np.loadtxt(path, skiprows=0, usecols=(0, 1))  # Petra/DNA_Paint

            if plot:
                fig, ax = plt.subplots()
                ax.plot(self.points[:, 0], self.points[:, 1], 'k.', markersize=0.8, alpha=0.5); ax.hold(True)
                ax.set_xlim(self.points[:, 0].min(), self.points[:, 0].max())
                ax.set_ylim(self.points[:, 1].min(), self.points[:, 1].max())
                ax.set_aspect('equal', adjustable='box'); ax.hold(False)


def saveparameters(file_name, dict1={}, dict2={}, dict3={}):

    f = open(file_name, 'wb')
    for ii in range(len(dict1)):
        name = dict1.keys()[ii]
        value = str(dict1[name])
        f.write("%s \t %s\n" % (name, value))

    for ii in range(len(dict2)):
        name = dict2.keys()[ii]
        value = str(dict2[name])
        f.write("%s \t %s\n" % (name, value))
    for ii in range(len(dict3)):
        name = dict3.keys()[ii]
        value = str(dict3[name])
        f.write("%s \t %s\n" % (name, value))

    f.close()


def generate_pattern_rnd(kwargs, num=1, save_data={}):
    """
    This function generates a synthetic point pattern (single-particle synthetic data).

    Input:
    -----------------
    cell_size (float): cell size in mic2
    distribution (string): probability distribution of the clusters, either 'uniform' or 'gaussian'
    geometry (string): geometrical shape of the clusters, either 'circle' or ...
    mean_radius (float): mean value of the radius normal pdf
    std_radius (float): std of the radius normal pdf
    density_inside (float): points/nm2 [see Le15]
    enrichment_ration (float): denisty inside / denisty background

    Output:
    -----------------
    points (2d-array): point pattern, X[:, 0] (x-coordinate) and X[_, 1] (y-coordinate). Returns the last pp of num.

    """
    import sys

    cell_size = kwargs.get('cell_size', 10*10**6)
    cell_width = np.sqrt(cell_size)
    n_points_total = kwargs.get('n_points_total', None)
    n_clusters = kwargs.get('n_clusters', 1)
    distribution = kwargs.get('distribution', 'uniform')
    geometry = kwargs.get('geometry', 'circle')
    mean_radius = kwargs.get('mean_radius', 50)
    std_radius = kwargs.get('std_radius', 0)
    enrichment_ratio = kwargs.get('enrichment_ratio', 20)

    for ii in range(0, num):
        print 'Generate pattern %d...' % ii
        area_c = np.pi*mean_radius**2
        area_b = (cell_size-n_clusters*area_c)

        n_cb = enrichment_ratio*area_c/area_b  # N_{cluster}/N_{background}
        n_points_background = int(1/(1+n_clusters*n_cb)*n_points_total)
        n_points_clusters = int((n_points_total-n_points_background) / n_clusters)  # same number of points/cluster
        n_points_in_cluster = int(n_clusters*n_points_clusters)
        n_points_total = n_points_background + n_points_in_cluster

        density_b = n_points_background/area_b
        density_c = enrichment_ratio*density_b

        print 'radius cluster = %.2f ' % mean_radius
        print 'number of points in one cluster = %.2f number of points in background' % n_cb
        print 'density background = %.0f points/u2' % density_b
        print 'density cluster = %.0f points/u2' % density_c
        print 'R = inside density/density outside = %.0f ' % enrichment_ratio

        print 'num. points total = %d points' % n_points_total
        print 'num. points in background = %d points' % n_points_background
        print 'num. points in (all) clusters = %d points' % n_points_in_cluster
        print 'num. points in (one) cluster (mean) = %d points' % np.mean(n_points_clusters)

        radius = np.vstack([np.matlib.repmat(ii, n_points_clusters, 1)
                           for ii in np.random.normal(loc=mean_radius, scale=std_radius,
                                                      size=(n_clusters,))]).reshape((n_points_in_cluster,))

        if distribution is 'uniform':  # same as in (Levet et al., 2015).
            points_background = np.random.uniform(low=0, high=cell_width, size=(n_points_background, 2))

            rnd_theta = np.random.uniform(low=0, high=(1.-sys.float_info.epsilon)*2*np.pi, size=(n_points_in_cluster,))
            rnd_radius = np.multiply(radius, np.sqrt(np.random.uniform(low=0, high=1, size=(n_points_in_cluster,))))
            clusters = np.zeros(shape=(n_points_in_cluster, 2))
            clusters[..., 0] = np.multiply(rnd_radius, np.cos(rnd_theta))
            clusters[..., 1] = np.multiply(rnd_radius, np.sin(rnd_theta))

    #        centers_clusters = np.random.uniform(low=0+np.min(radius)+3*std_radius, high=cell_width-np.max(
    #                                             radius)-3*std_radius, size=(n_clusters, 2))
            centers_clusters = [[cell_width*0.5, cell_width*0.5]]
            rep_centers = np.vstack([np.matlib.repmat(centers_clusters[ii], n_points_clusters, 1)
                                     for ii in range(len(centers_clusters))])
            clusters = rep_centers + clusters

            points = np.zeros(shape=(n_points_total+1, 2))
            points[0:n_points_background] = points_background  # np.zeros(shape=points_background.shape)
            points[n_points_background:-1] = clusters

            plt.figure(); plt.plot(points[..., 0], points[..., 1], 'k.', markersize=0.5)
            plt.show(); plt.axis('equal')

        if distribution is 'normal':  # refers to gaussian clusters as in (Delanchy et al., 2015)
            from sklearn.datasets.samples_generator import make_blobs

            X, y_true = make_blobs(n_samples=400, centers=4, cluster_std=0.60, random_state=0); X = X[:, ::-1]
            plt.scatter(X[:, 0], X[:, 1])

            from sklearn.cluster import KMeans

            kmeans = KMeans(4, random_state=0)
            labels = kmeans.fit(X).predict(X)
            plt.figure(); plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')

            xx = np.array([-0.51, 51.2])
            yy = np.array([0.33, 51.6])
            means = [xx.mean(), yy.mean()]
            stds = [xx.std() / 3, yy.std() / 3]
            corr = 0.8         # correlation
            covs = [[stds[0]**2, stds[0]*stds[1]*corr],
                    [stds[0]*stds[1]*corr,           stds[1]**2]]

            m = np.random.multivariate_normal(means, covs, 1000).T
            plt.scatter(m[0], m[1])

        print(save_data.get('output_dir') + r'syntheticpp_cw%d_enrichmentr%d_r%.2f_Nt%d_Ch1_%d.txt'
                       % (cell_size, enrichment_ratio, np.mean(n_points_clusters), n_points_total, ii))
        if save_data.get('save', 0):
            np.savetxt(save_data.get('output_dir') + r'syntheticpp_cw%d_enrichmentr%d_r%d_Nt%d_Ch1_%d.txt'
                       % (cell_size, enrichment_ratio, np.mean(n_points_clusters), n_points_total, ii), points)

    return points/float(kwargs.get('unit_size', 1.))


def generate_circularpattern_det(kwargs):
    """
    This function generates a synthetic point pattern (single-particle synthetic data).

    Input:
    -----------------
    cell_size (float): cell size in mic2
    distribution (string): probability distribution of the clusters, either 'uniform' or 'gaussian'
    geometry (string): geometrical shape of the clusters, either 'circle' or ...
    density_inside (float): points/nm2 [see Le15]
    enrichment_ration (float): denisty inside / denisty background

    Output:
    -----------------
    points (2d-array): point pattern, X[:, 0] (x-coordinate) and X[_, 1] (y-coordinate). Unit: STORM pixel size (
    1px=160nm)

    """

    import sys

    cell_size = kwargs.get('cell_size', 10*10**6)
    radius = kwargs.get('radius')
    n_clusters = radius.shape[0]
    distribution = kwargs.get('distribution', 'uniform')
    density_inside = kwargs.get('density_inside', 0.01)
    enrichment_ratio = kwargs.get('enrichment_ratio', 20)

    cell_width = np.sqrt(cell_size)

    n_points_background = int(density_inside/enrichment_ratio*cell_size)

    temp = 1./enrichment_ratio*100
    n_points_clusters = (density_inside*np.pi*radius**2).astype(int)
    n_points_no_cluster = np.sum(n_points_clusters)
    n_points_total = n_points_background + n_points_no_cluster

    radius = np.vstack([np.matlib.repmat(radius[ii], n_points_clusters[ii], 1)
                        for ii in range(len(radius))]).reshape((n_points_no_cluster,))

    print 'outside density = %.0f [percent] inside density' % temp
    print 'n_points_background = %.3f points' % n_points_background
    print 'n_point_total = %d points' % n_points_total
    print 'n_points_clusters (mean) = %d points/cluster' % np.mean(n_points_clusters)

    if distribution is 'uniform':  # same as in (Levet et al., 2015).
        points_background = np.random.uniform(low=0, high=cell_width, size=(n_points_background, 2))

        rnd_theta = np.random.uniform(low=0, high=(1.-sys.float_info.epsilon)*2*np.pi, size=(n_points_no_cluster,))
        rnd_radius = np.multiply(radius, np.sqrt(np.random.uniform(low=0, high=1, size=(n_points_no_cluster,))))
        clusters = np.zeros(shape=(n_points_no_cluster, 2))
        clusters[..., 0] = np.multiply(rnd_radius, np.cos(rnd_theta))
        clusters[..., 1] = np.multiply(rnd_radius, np.sin(rnd_theta))

        centers_clusters = np.random.uniform(low=0+np.min(radius), high=cell_width-np.max(
                                             radius), size=(n_clusters, 2))
        rep_centers = np.vstack([np.matlib.repmat(centers_clusters[ii], n_points_clusters[ii], 1) for ii in range(len(
                                centers_clusters))])
        clusters = rep_centers + clusters

        points = np.zeros(shape=(n_points_total+1, 2))
        points[0:n_points_background] = points_background  # np.zeros(shape=points_background.shape)
        points[n_points_background:-1] = clusters

        plt.figure(); plt.plot(points[..., 0], points[..., 1], 'k.', markersize=0.5)
        plt.show(); plt.axis('equal')

    if distribution is 'normal':  # refers to gaussian clusters as in (Delanchy et al., 2015)
        from sklearn.datasets.samples_generator import make_blobs

        X, y_true = make_blobs(n_samples=400, centers=4, cluster_std=0.60, random_state=0); X = X[:, ::-1]
        plt.scatter(X[:, 0], X[:, 1])

        from sklearn.cluster import KMeans

        kmeans = KMeans(4, random_state=0)
        labels = kmeans.fit(X).predict(X)
        plt.figure(); plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')

        xx = np.array([-0.51, 51.2])
        yy = np.array([0.33, 51.6])
        means = [xx.mean(), yy.mean()]
        stds = [xx.std() / 3, yy.std() / 3]
        corr = 0.8         # correlation
        covs = [[stds[0]**2, stds[0]*stds[1]*corr],
                [stds[0]*stds[1]*corr,           stds[1]**2]]

        m = np.random.multivariate_normal(means, covs, 1000).T
        plt.scatter(m[0], m[1])

    return points/160.


def generate_image():

    """
    This function generates a synthetic image.

    """
    radius = 0.5*np.array([5, 14])
    # radius = 0.5*np.array([28, 33, 39])
    # radius = 0.5*np.arange(11, 23)
    radius_n = np.array([1, 1])
    # radius_n = np.array([1, 3])
    # radius_n = np.arange(1, len(radius)+1)
    maxradius = np.max(radius)
    width = int(2*2*maxradius*np.max(radius_n))
    height = int(2*1.9*maxradius*radius.shape[0])

    print 'diameter:'
    print 2*radius
    print 'number of blobs = %d' % np.sum(radius_n)
    # color = np.linspace(100, 255, num=np.max(radius_n), dtype=float)
    image = np.zeros(shape=(height, width))

    xv, yv = np.meshgrid(np.linspace(0, width, width+1), np.linspace(0, height, height+1))
    print xv

    grid_centers_y = maxradius + np.arange(maxradius, height-2*maxradius, 3*maxradius)
    grid_centers_x = maxradius + np.arange(maxradius, width-2*maxradius, 3*maxradius)

    for ii, itemii in enumerate(radius):
        for jj in range(radius_n[ii]):
            pos = np.where(((xv-grid_centers_x[jj])**2 + (yv-grid_centers_y[ii])**2) < radius[ii]**2)
            image[pos] = 255  # color[jj]

    return image


def violin_plot(ax, data1, pos1, data2={}, pos2={}, bp=False, xlabel='', ylabel=''):
    """
    create violin plots on an axis
    """

    from scipy.stats import gaussian_kde

    pos = np.arange(0, len(pos1))
    dist = max(pos)-min(pos)

    w = min(0.15*max(dist, 1.0), 0.5)
    print data1
    for p, d in enumerate(data1):
        k = gaussian_kde(d)  # calculates the kernel density
        m = k.dataset.min()  # lower bound of violin
        M = k.dataset.max()  # upper bound of violin
        x = np.arange(m, M, (M-m)/100.)  # support for violin
        v = k.evaluate(x)  # violin profile (density curve)
        v = v/v.max()*w  # scaling the violin to the available space
        ax.fill_betweenx(x, p, v+p, facecolor='lightgray', edgecolor='gray', alpha=0.3)
        ax.fill_betweenx(x, p, -v+p, facecolor='lightgray', edgecolor='gray', alpha=0.3)

    # dist = max(pos2) - min(pos2)
    # w = min(0.15 * max(dist, 1.0), 0.5)
    # for d, p in zip(data2, pos2):
    #     k = gaussian_kde(d)  # calculates the kernel density
    #     m = k.dataset.min()  # lower bound of violin
    #     M = k.dataset.max()  # upper bound of violin
    #     x = np.arange(m, M, (M - m) / 100.)  # support for violin
    #     v = k.evaluate(x)  # violin profile (density curve)
    #     v = v / v.max() * w  # scaling the violin to the available space
    #     ax.fill_betweenx(x, p, v + p, facecolor='sage', edgecolor='gray', alpha=0.3)
    #     ax.fill_betweenx(x, p, -v + p, facecolor='sage', edgecolor='gray', alpha=0.3)

    if bp:
        # ax.boxplot(data1.transpose(),notch=1,positions=pos1,vert=1)
        # plt.hold(True)
        # ax.boxplot(data2.transpose(), notch=1, positions=[1,3,5], vert=1)
        plt.xticks(pos, pos1)
        # plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], size=16)
        # plt.xlim([-0.5,5.5])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)


def plot_frameno(points, frame_no):
    import matplotlib.colors as colors
    import matplotlib.cm as cmx

    cmap = plt.cm.jet
    cNorm = colors.Normalize(vmin=np.min(frame_no), vmax=np.max(frame_no))  # vmax=values[-1]) . LogNorm,
    # Normalize
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)  # norm=cNorm
    scalarMap.set_array(frame_no)

    fig, ax = plt.subplots()

    for ii in range(frame_no.shape[0]):
        deviation = frame_no[ii]
        blob_color = scalarMap.to_rgba(deviation)
        ax.plot(points[ii][0], points[ii][1], marker='o', fillstyle='full', color=blob_color, markersize=3,
                markeredgewidth=0.0)
        plt.hold(True)
    plt.colorbar(scalarMap, label='Frame number')
    plt.hold(False)


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    import matplotlib as mpl
    if N < 8:
        colormap_rgb = [(0.9,0.62,0), (0.34,0.7,0.91), (0, 0.62, 0.45), (0.94,0.89,0.25), (0,0.45,0.7),
                                          (0.83,0.368,0), (0.8,0.474,0.654)]
        colormap = mpl.colors.ListedColormap(colormap_rgb[0:N])
    else:
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm

        colormap = cm.viridis
        normalize = mcolors.Normalize(vmin=0, vmax=N-1)
        color = colormap(normalize(range(N)))
        colormap = mpl.colors.ListedColormap(color)
    return colormap


def permute_labels(siftclusters, k, centers_permuted0):
    """
    This function assignes a permuted version of the labels (or classes) just to be consistent with the colorcode of
    the previous iteration in k (number of clusters, words,...)
    """

    labels_permuted = np.full_like(siftclusters.labels_, -1)
    centers_permuted = np.full_like(siftclusters.cluster_centers_, -1)
    list_labels = []
    for l0, c in enumerate(centers_permuted0):
        l1 = np.argmin(np.linalg.norm(c - siftclusters.cluster_centers_, axis=1) ** 2)
        labels_permuted[np.where(siftclusters.labels_ == l1)] = l0
        centers_permuted[l0] = siftclusters.cluster_centers_[l1]
        list_labels.append(l1)
    labels_permuted[np.where(labels_permuted == -1)] = k - 1
    centers_permuted[k - 1] = siftclusters.cluster_centers_[np.setdiff1d(range(k), list_labels)[0]]
    siftclusters.labels_ = labels_permuted
    siftclusters.cluster_centers_ = centers_permuted
    centers_permuted0 = centers_permuted

    return centers_permuted0, siftclusters


def select_labels_image(feature_all, image_no):

    feature = feature_all[image_no]
    siftclusters_labels_image = [feature['word_label'][ii][jj] for ii, orientation in
                                enumerate(feature['orientation']) for jj, ori in enumerate(orientation)]

    return siftclusters_labels_image


def points_2dcrop(points, rangexy):
    """
    select smaller region

    Input:
    -------------
    points (nx2 array)
    rangexy (list of floats): [x0, x1, y0, y1] cropped area is [x0, x1]x[y0,y1]

    Output:
    --------------
    new_points (nx2 array)
    """

    x0 = rangexy[0]; x1 = rangexy[1]; y0 = rangexy[2]; y1 = rangexy[3]
    temp_pos = np.where(points[:, 0] > x0)
    new_points = points[temp_pos]
    temp_pos = np.where(new_points[:, 0] < x1)
    new_points = new_points[temp_pos]

    temp_pos = np.where(new_points[:, 1] > y0)
    new_points = new_points[temp_pos]
    temp_pos = np.where(new_points[:, 1] < y1)
    new_points = new_points[temp_pos]

    return new_points


def remove_points_rnd(points, percent=10):

    """
    Randomly remove points by percent [%]
    """
    num_points = points.shape[0]
    rnd_pos = np.random.randint(num_points-1, size=int(num_points*(1-1/100.*percent)))

    new_points = points[rnd_pos]

    return new_points


def memory_usage():

    import os

    mem = str(os.popen('free -t -m').readlines())
    """
    Get a whole line of memory output, it will be something like below
    ['             total       used       free     shared    buffers     cached\n',
    'Mem:           925        591        334         14         30        355\n',
    '-/+ buffers/cache:        205        719\n',
    'Swap:           99          0         99\n',
    'Total:        1025        591        434\n']
     So, we need total memory, usage and free memory.
     We should find the index of capital T which is unique at this string
    """
    T_ind = mem.index('T')
    """
    Than, we can recreate the string with this information. After T we have,
    "Total:        " which has 14 characters, so we can start from index of T +14
    and last 4 characters are also not necessary.
    We can create a new sub-string using this information
    """
    mem_G = mem[T_ind + 14:-4]
    """
    The result will be like
    1025        603        422
    we need to find first index of the first space, and we can start our substring
    from from 0 to this index number, this will give us the string of total memory
    """
    S1_ind = mem_G.index(' ')
    mem_T = mem_G[0:S1_ind]
    """
    Similarly we will create a new sub-string, which will start at the second value.
    The resulting string will be like
    603        422
    Again, we should find the index of first space and than the
    take the Used Memory and Free memory.
    """
    mem_G1 = mem_G[S1_ind + 8:]
    S2_ind = mem_G1.index(' ')
    mem_U = mem_G1[0:S2_ind]

    mem_F = mem_G1[S2_ind + 8:]
    print 'Summary = ' + mem_G
    print 'Total Memory = ' + mem_T + ' MB'
    print 'Used Memory = ' + mem_U + ' MB'
    print 'Free Memory = ' + mem_F + ' MB'

    CPU_Pct = str(round(float(
        os.popen('''grep 'cpu ' /proc/stat | awk '{usage=($2+$4)*100/($2+$4+$5)} END {print usage }' ''').readline()),
                        2))

    # print results
    print("CPU Usage = " + CPU_Pct)


def read_features(input_files, redo_trainingset=0, training_size=0.7):

    import pickle

    feature_all_aux = []; file_dirs = []; trainingset_filenames = []; trainingset_filenames_all = []
    for kk in range(len(input_files)):
        print '\t reading output variables in', '.../' + '/'.join(input_files[kk].split('/')[-3:])
        file_variables = open(input_files[kk],'rb')
        output = pickle.load(file_variables); file_variables.close()
        feature_all = output[0]; file_dir = output[1]; dict_inputfile = output[2]
        dict_image = output[3]; dict_sift = output[4]

        feature_all_aux.append(feature_all)

        if len(output) < 6:
            print '\t\tThis dataset does not contain a training set. Generating training set with', training_size*100, '% of dataset...'
            trainingset_filenames = generate_training_set(feature_all, size_set=training_size)
        else:
            print '\t\tThis dataset contains already a training set.'
            trainingset_filenames = output[5]
            if redo_trainingset:
                print '\t\tOverwriting... generating training set with', training_size * 100, '% of dataset...'
                trainingset_filenames = generate_training_set(feature_all, size_set=training_size)
        if len(output) < 6 or redo_trainingset:
            print '\t\tNew trainingset text file saved in variables.'
            file_variables = open(input_files[kk], 'wb')
            pickle.dump([output[0], output[1], output[2], output[3], output[4], trainingset_filenames],
                        file_variables)
            file_variables.close()

        output = []; del output
        trainingset_filenames_all.extend(trainingset_filenames)
        file_dirs.append(file_dir)

    feature_all = []; del feature_all
    feature_all = [feature for exp in feature_all_aux for feature in exp]; feature_all_aux = []; del feature_all_aux

    return feature_all, file_dirs, dict_inputfile, dict_image, dict_sift, trainingset_filenames_all


def train_test_split(feature_all, trainingset_filenames):

    feature_all_test = []
    if len(trainingset_filenames) > 0:
        feature_all_aux = []; feature_all_test_aux = []
        print '\t\t\t- size of training set = ', len(trainingset_filenames), 'of', len(feature_all)
        print '\t\tTraining set:'
        for ii, feat in enumerate([feature for feature in feature_all if feature['file_name'] in
                trainingset_filenames]):
                print '\t\t\t', ii, ''.join(feat['file_dir'].split('/')[-2:]), '\t', feat['file_name']
        print '\t\tTest set:'
        for ii, feat in enumerate([feature for feature in feature_all if feature['file_name'] not in \
                trainingset_filenames]):
                print '\t\t\t', ii, ''.join(feat['file_dir'].split('/')[-2:]), '\t', feat['file_name']
        # for file_name in [feature['file_name'] for feature in feature_all if feature['file_name'] \
        #         not in trainingset_filenames]:
        #     print '\t\t\t\t', file_name

        feature_all_aux.append(
            [feature for feature in feature_all if feature['file_name'] in trainingset_filenames])
        feature_all_test_aux.append(
            [feature for feature in feature_all if feature['file_name'] not in trainingset_filenames])

        feature_all = []; del feature_all
        feature_all = [feature for exp in feature_all_aux for feature in exp];
        feature_all_aux = []; del feature_all_aux

        feature_all_test = [feature for exp in feature_all_test_aux for feature in exp]
        feature_all_test_aux = []; del feature_all_test_aux
    else:
        print '\tno training/test split is performed.'

    return feature_all, feature_all_test


def save_features_small(input_files):

    import pickle
    for kk in range(len(input_files)):
        print '\t reading output variables in', input_files[kk]
        file_variables = open(input_files[kk], 'rb')
        output = pickle.load(file_variables)
        feature_all = output[0]; file_dir = output[1]; dict_inputfile = output[2]
        dict_image = output[3]; dict_sift = output[4]
        if len(output) < 6:
            trainingset_filenames = generate_training_set(feature_all)
        else: trainingset_filenames = output[5]
        file_variables.close()
        output = []; del output
        print '\t removing vor and image and saving to', input_files[kk]+'_small'
        for feature in feature_all:
            print '\t\t\t', feature['file_name']
            feature.pop('vor', None); feature.pop('image', None)
        file_variables = open(input_files[kk]+'_small', 'wb')
        pickle.dump([feature_all, file_dir, dict_inputfile, dict_image, dict_sift, trainingset_filenames],
                    file_variables)
        file_variables.close()

    return 1


def generate_training_set(feature_all, size_set=0.8):

    import random
    trainingset = random.sample(range(len(feature_all)), k=int(size_set*len(feature_all)))
    trainingset.sort()
    trainingset_filenames = [feature['file_name'] for feat, feature in enumerate(feature_all) if feat in
                                 trainingset]
    return trainingset_filenames


def colorbar(mappable):

    """
    function that resized colorbar to hight of the figure. Taken from https://joseph-long.com/writing/colorbars/

    call: img = ax.imshow(data); colorbar(img)

    """

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def make_meshgrid(x, h=.02):
    """Create a mesh of points to plot in

    Parameters
    ----------
    x: array with n_sample \times n_features dimension to base meshgrid on [x0, x1, x2, ...].
    h: stepsize for meshgrid, optional

    Returns
    -------
    ndarrays
    """

    mesh = []
    h = 1./100*abs(x[:, 0].max()+abs(0.5*x[:, 0].max())-x[:, 0].min()-abs(0.5*x[:, 0].min()))
    print h

    for ii in range(x.shape[1]):  # in dimension - each coordinate
        mesh.append(np.arange(x[:, ii].min()-abs(0.5*x[:, ii].min()), x[:, ii].max()+abs(0.5*x[:, ii].max()), h))

    return np.meshgrid(*mesh)

    # x_min, x_max = x.min() - 1, x.max() + 1
    # y_min, y_max = y.min() - 1, y.max() + 1
    # xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
    #                      np.arange(y_min, y_max, h))
    # return xx, yy


def plot_contours(ax, svm_model, xx, yy, **params):
    """Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """

    Z = svm_model.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)

    return out


def find_file_dirs(feature_all):

    if not isinstance(feature_all, list): feature_all = [feature_all]  # only one sample

    file_dirs_all = np.asarray([feature['file_dir'] for feature in feature_all])
    file_dirs = file_dirs_all[np.sort(np.unique(file_dirs_all, return_index=True)[1])]

    return file_dirs
