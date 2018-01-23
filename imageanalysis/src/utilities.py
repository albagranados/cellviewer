import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
# configure latex plots
matplotlib.rcParams['text.usetex'] = True  # uselatex labels
matplotlib.rcParams['font.family'] = 'serif'


def read_points(fileDir, fileName=None, fileExt='.txt', storm=1, channels_num=2, out_channel='all',
                save_output_dir=None, plot=1):
    """
    This function reads point pattern(s) stored in fileDir

    Input:
    ----------------------------------
    fileDir (string): file(s) directory
    fileName (string): if not None, then read multple data sets in fileDir
    fileExt (string): extension of the pointpattern files
    storm (boolean): if the datafile is in storm format (precision, coordinates...)
    save_output_dir (string): if not None, then save the point pattern(s) in a txt files as a nx2 matrix in the output
                                directory (string)

    Output:
    ----------------------------------
    points: points (if multiple data sets, the last one) stored in a nx2 matrix

    """

    if fileName is not None:
        fileNames = [fileName]
    else:
        fileNames = os.listdir(fileDir)

    for ii in range(len(fileNames)):
        fileName = fileNames[ii].split(fileExt)[0]  # '2017-12-14_HeLa_DMSO_000_list_m1000DC_CTS_0851-1731_allChs'
        print('Reading data set: ', fileName)

        path = os.path.join(fileDir, fileName + fileExt)  # if we want to run in the Python Console

        if storm:
            points = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=(3, 4))  # export x_c and y_c
            channels = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=(0,))  # export x_c and y_c
            # area = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=6)  # export x_c and y_c
            # width = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=7)  # export x_c and y_c
            # aspect_ratio = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=9)
            # bg = np.loadtxt(path, skiprows=1, delimiter='\t', usecols=10)

        if not storm:
            points = np.loadtxt(path, skiprows=0, usecols=(0, 1))  # Petra/DNA_Paint

        if save_output_dir is not None:
            if channels_num == 2:
                np.savetxt(save_output_dir + fileName + '_Ch1' + fileExt, points[np.where(channels == 1)])
                np.savetxt(save_output_dir + fileName + '_Ch2' + fileExt, points[np.where(channels == 2)])
            if channels_num == 1:
                np.savetxt(save_output_dir + fileName + fileExt, points)

        if out_channel is not 'all':
            points = points[np.where(channels == out_channel)]
            print 'output points in channel %d' %out_channel

        if plot:
            fig, ax = plt.subplots()
            ax.plot(points[:, 0], points[:, 1], 'k.', markersize=2); plt.show(); plt.hold(True)

            ax.set_xlim(points[:, 0].min(), points[:, 0].max()); ax.set_ylim(points[:, 1].min(), points[:, 1].max())
            ax.set_aspect('equal', adjustable='box')

            plt.hold(False)

    return points


def saveparameters(fileName, dict1={}, dict2={}, dict3={}):

    f = open(fileName, 'wb')
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
        print 'Generate pattern %d...' %ii
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
            covs = [[stds[0]**2 , stds[0]*stds[1]*corr],
                    [stds[0]*stds[1]*corr,           stds[1]**2]]

            m = np.random.multivariate_normal(means, covs, 1000).T
            plt.scatter(m[0], m[1])

        print(save_data.get('output_dir') + r'syntheticpp_cw%d_enrichmentr%d_r%.2f_Nt%d_Ch1_%d.txt'
                       % (cell_size, enrichment_ratio, np.mean(n_points_clusters), n_points_total, ii))
        if save_data.get('save',0):
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
        covs = [[stds[0]**2 , stds[0]*stds[1]*corr],
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
    print 'number of blobs = %d' %np.sum(radius_n)
    color = np.linspace(100, 255, num=np.max(radius_n), dtype=float)
    image = np.zeros(shape=(height, width))

    xv, yv = np.meshgrid(np.linspace(0, width, width+1), np.linspace(0, height, height+1))
    print xv

    grid_centers_y = maxradius + np.arange(maxradius, height-2*maxradius, 3*maxradius)
    grid_centers_x = maxradius + np.arange(maxradius, width-2*maxradius, 3*maxradius)

    for ii, itemii in enumerate(radius):
        for jj in range(radius_n[ii]):
            pos = np.where(((xv-grid_centers_x[jj])**2 + (yv-grid_centers_y[ii])**2) < radius[ii]**2)
            image[pos] = 255  #color[jj]

    return image
