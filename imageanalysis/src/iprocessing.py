import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams["text.usetex"] = True; matplotlib.rcParams['font.family'] = 'serif'  # configure latex plots


def pattern2image(points, pixel_size):
    """
    Converts a point pattern (spatial object) to a pixel image. The value in each pixel is the number of points
    falling in that pixel. The function counts the number of points falling in each pixel. Kinda "quadrat analysis"
    in GIS.
    The pixel size, in nm, to use to perform the cluster analysis. A 2D histogram of the number of localizations per
    pixel could be created

    Input:
    ---------
    points (ndarray): input spatial point pattern. Is equivalent to scipy.spatial.Voronoi points

    Output:
    ---------
    image (2d-array): is the row_resolution image corresponding to the given point pattern
    image_ptslabel (2d-array): each (row,col) contains a 1d-array of the labels of the points within that pixel.
    Used at image2pattern for a given (roi) image. empty array [] corresponds to pixel with no points
    """

    import math

    points_x = points[:, 0]  # cartesian x-coordinate
    points_y = points[:, 1]  # cartesian y-coordinate
    num_points = points_x.size  # num. of points in the original spatial point pattern

    range_x = points_x.ptp()  # Peak to peak (maximum - minimum) value along a given axis.
    range_y = points_y.ptp()

    nx = int(math.ceil(range_x / pixel_size))  # num. of pixels in the x direction
    ny = int(math.ceil(range_y/pixel_size))  # num. of pixels in the y direction

    min_pointx = np.min(points_x)
    min_pointy = np.min(points_y)
    grid_x = np.linspace(min_pointx, min_pointx + nx * pixel_size, nx + 1)
    grid_y = np.linspace(min_pointy, min_pointy + ny * pixel_size, ny + 1)

    image = np.zeros((nx, ny))  # np.zeros( (nx+1, ny+1) )
    image_ptslabel = np.ndarray(shape=image.shape, dtype=object)

    count = 0
    pos = np.linspace(0, num_points - 1, num=num_points, dtype=int)  # global points labelling
    pos_out = np.copy(pos)  # global
    for ii, pixel_x in enumerate(grid_x[1:]):
        pos_in_x = np.where(points_x[pos_out] <= pixel_x)  # global
        pos_in = pos_out[pos_in_x]
        # print 'X ----',  pixel_x,  "  ", pos_in
        # if any(x<=0 for x in points_x[pos_out]-pixel_x): print points_x[pos_in]
        pos_out = np.delete(pos_out, pos_in_x)
        for jj, pixel_y in enumerate(grid_y[1:]):
            # print "Y ***", pixel_y, "  ", points_y[pos_in]
            pos_in_xy = np.where(points_y[pos_in] <= pixel_y)  # local in pos_in
            # if any(x <= 0 for x in points_y[pos_in] - pixel_y):  print points_y[pos_in][pos_in_xy]
            image[ii][jj] = len(pos_in_xy[0])
            image_ptslabel[ii][jj] = pos_in[pos_in_xy[0]]
            count += len(pos_in_xy[0])
            # print "count: ", image[j][i]
            pos_in = np.delete(pos_in, pos_in_xy)

    return image, image_ptslabel


def image2pattern(image_roi, points, image_ptslabel):
    """
    Given a (roi) low-resolution image, the function returns to the original point pattern within that roi

    Input:
    ---------
    image_roi (2d-array): binary image. The (roi) low-resolution image pattern
    points (nd-array): full original spatial point pattern. Is equivalent to scipy.spatial.Voronoi points
    image_ptslabel (2d-array): each (row,col) contains a 1d-array of the labels of the points within that pixel.

    Output:
    ---------
    roi_pattern (2d-array): is a point pattern containing the points within the (roi) input low-reso image
    """

    points_x = points[:, 0]  # cartesian x-coordinate
    points_y = points[:, 1]  # cartesian y-coordinate

    pixels_roi = np.where(image_roi)
    aux = [xx.tolist() for xx in image_ptslabel[pixels_roi] if xx.tolist() != []]
    roi_ptslabel = [label for pixels in aux for label in pixels]

    points_roi = np.ndarray(shape=(len(roi_ptslabel), 2), dtype=float)
    points_roi[:, 0] = points_x[roi_ptslabel]
    points_roi[:, 1] = points_y[roi_ptslabel]

    return points_roi


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


def image_crop(image, x0, x1, y0, y1):
    """
    select smaller region
    """
    new_image = image[x0:x1, y0:y1]

    return new_image


def plot_image(image, cmap='gray', interpolation='none', norm=None, plot_axis='on', hold=False):
    """
    Plot the image obtained with pattern2image... and more! NB: It plots the transpose as image[x][y]=points[x][y] (
    x-row; y-col)

    Input:
    ---------
    points_image: image obtained from the the tranformation pattern to image

    """
    from matplotlib.colors import LogNorm

    fig, ax = plt.subplots()
    if norm is 'log':
        plt.imshow(image.T, interpolation=interpolation, cmap=cmap, origin='lower', norm=LogNorm())
        # plt.colorbar()
    else:
        plt.imshow(image.T, interpolation=interpolation, cmap=cmap, origin='lower')

    # plt.colorbar()
    fig.hold(1)
    if plot_axis is 'off':
        plt.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    fig.hold(hold)
    return fig, ax


def blur_image(image, t):
    """
    Blurs image with a Gaussian filter with variance t

    Input:
    -------
    t (float): variance of the Gaussian kernel

    Output:
    ------
    image_blurred (ndarray): image with applied filter
    """

    from scipy import ndimage

    gauss_kernels = get_gauss_filter(t)
    g = gauss_kernels.get('g')

    image_blurred = ndimage.convolve(image, g[None, ...], mode='constant', cval=0.0)
    image_blurred = ndimage.convolve(image_blurred, g[..., None], mode='constant', cval=0.0)

    return image_blurred


def find_feature(image, kwargs, image_descr=None):
    """
    Find ROI corresponding to the nucleus to be processed.

    This function requires pattern2image.

    Input:
    ---------
    image (ndarray)
    kwargs (dictionary dict_roi) with:
        't' (float): the variance of the Gaussian kernels.
        'feature_name' (string): 'edge' or 'blob'
        'thresholding' (bool) - True or False, depending whether we want to capture only the strongest feature responses
        'threshold_percent' (float 0<.<=1): the threshold_percent*100 strongest responses. Default 0.4.
    image: if not empty, then we use image_descr for the feature descriptors, different from the feature detection.

    Output:
    ---------
    edgemat (dictionary)

    """
    feature_name = kwargs.get('feature_name', 'blob')
    # t = kwargs.get('t', 10)
    # thresholding = kwargs.get('thresholding', False)
    # nscales = kwargs.get('nscales', 1)

    if feature_name == 'edge':
        feature = get_edge(image, kwargs)

    if feature_name == 'blob':
        feature = get_blob(image, kwargs, image_descr=image_descr)

    #  # begin: this code below has been moved to get_blob -> we compute orientations just in the thresholded features.
    # argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    # strength = feature.get('strength')  # 2d-array
    # scale = feature.get('scale')
    # if thresholding:
    #     threshold_percent = kwargs.get('threshold_percent', 0.4)  # default 40% down the maximum response
    #     threshold = np.max(strength[argmaxgrad])-threshold_percent*(np.max(strength[argmaxgrad]) -
    #                                                                        np.min(strength[argmaxgrad]))
    # else:
    #     threshold = -float('inf')
    #
    # feature['threshold'] = threshold
    #
    # temp = np.where(strength[argmaxgrad] > threshold)  # tuple
    # argmaxgrad_threshold = (argmaxgrad[0][temp[0]], argmaxgrad[1][temp[0]])
    # feature['argmaxgrad'] = argmaxgrad_threshold
    # scale_threshold = scale[temp[0]]
    # feature['scale'] = scale_threshold
    # image_roi = np.zeros(image.shape)
    # image_roi[argmaxgrad_threshold] = 1
    # image_roi = flood_fill(image_roi)
    # feature['image_roi'] = image_roi

    # orientation = feature.get('orientation', [])
    # histogram_descr = feature.get('histogram_descr', [])
    # if len(orientation) > 0:
    #     feature['orientation'] = [ori for ii, ori in enumerate(orientation) if ii in temp[0]]  # orientation[temp[0]]
    # if len(histogram_descr) > 0:
    #     feature['histogram_descr'] = [hist for ii, hist in enumerate(histogram_descr) if ii in temp[0]]
    # # end

    return feature


def flood_fill(test_array, h_max=255):
    """
    Function equivalent to imfill in Matlab. Fills the closed surface.
    See http://stackoverflow.com/questions/36294025/python-equivalent-to-matlab-funciton-imfill-for-grayscale

    """
    from scipy import ndimage

    input_array = np.copy(test_array)
    el = ndimage.generate_binary_structure(2, 2).astype(np.int)
    inside_mask = ndimage.binary_erosion(~np.isnan(input_array), structure=el)
    output_array = np.copy(input_array)
    output_array[inside_mask] = h_max
    output_old_array = np.copy(input_array)
    output_old_array.fill(0)
    el = ndimage.generate_binary_structure(2, 1).astype(np.int)
    while not np.array_equal(output_old_array, output_array):
        output_old_array = np.copy(output_array)
        output_array = np.maximum(input_array, ndimage.grey_erosion(output_array, size=(3, 3), footprint=el))
    return output_array


def get_blob(image, kwargs, image_descr=None):
    """
    detection of very Large OBject

    Input:
    ---------
    image (ndarray)
    t (float): the variance of the Gaussian kernels.

    Output:
    ---------
    dictionary with
     'argmaxgrad' (tupla for 2d-array - image): argmaxgrad
     'strength' (2d-array): blobstrength
     'scale' (1d-array length argmaxgrad): selected scale parameter (variance of Gaussian kernel)

    """
    from scipy.ndimage.filters import maximum_filter

    t = kwargs.get('t', 1)
    original_pixel_size = kwargs.get('original_pixel_size', 1)
    scale_pixel_size = kwargs.get('scale_pixel_size', 1)
    pixel_size = scale_pixel_size * original_pixel_size  # = analysis pixel size
    nscales = kwargs.get('nscales', 1)
    scale_resolution = kwargs.get('scale_resolution', 1)  # only for odd spacing
    max_filter_width = kwargs.get('max_filter_width', 3)
    max_filter_depth = kwargs.get('max_filter_depth', None)
    scale_spacing = kwargs.get('scale_spacing', 'log')
    scale_ini = kwargs.get('scale_ini', 1)  # initial radius of search
    scale_end = kwargs.get('scale_end', 10)
    compute_orientation = kwargs.get('compute_orientation', False); orientation = []
    compute_sift_descr = kwargs.get('compute_sift_descr', False); histogram_descr = []
    plot_graphics = kwargs.get('plot_graphics', False)

    gamma = 1  # parameter normalized derivatives
    if nscales == 1:
        scalespace = compute_space_derivatives(image, t)
        lxx = scalespace.get('lxx')
        lyy = scalespace.get('lyy')
        laplacian = lxx + lyy
        strength = -t ** gamma * laplacian  # +t** -> find minima

        dila = np.ones(shape=(max_filter_width, max_filter_width), dtype=float)
        dila[max_filter_width/2, max_filter_width/2] = 0
        strength_dila = maximum_filter(strength, footprint=dila, mode='constant', cval=-float('inf'))
        argmaxgrad = np.where(strength > strength_dila)

        scale = np.ones(len(argmaxgrad[0])) * t
        scale_range = scale
    else:
        print '\tComputing scale range...'
        if scale_spacing is 'odd':
            if kwargs.get('scale_range_is') is 'nm':
                scale_ini = int(scale_ini/pixel_size)
                scale_end = int(scale_end/pixel_size)
                scale_resolution = np.ceil(scale_resolution/pixel_size)
                print '\t\tWarning: if scale_resolution [nm] < analysis_pixel_size [nm], then scale_resolution = ' \
                      'analysis_pixel_size:'
                print '\t\t\t\tscale_resolution = %.2f [nm] and analysis_pixel_size = %.2f [nm]' \
                      % (scale_resolution*pixel_size, pixel_size)
            scale_ini = np.ceil((scale_ini - 1) / 2.)
            scale_end = np.ceil((scale_end - 1) / 2.)
            if max_filter_depth is None: max_filter_depth = nscales + 1
            scale_range = (3**-2)*(2*np.arange(scale_ini, scale_end+1, step=scale_resolution) + 1) ** 2
            nscales = len(scale_range)
        elif scale_spacing is 'log': 
            print '\Warning in building scale range: with log spacing and arbitrary nscales, maximum in Laplacian ' \
                  'direction can give many false positive -> noisy! Better use odd- spacing with a resolution.'
            if kwargs.get('scale_range_is') is 'pixel':
                scale_range = (3**-2)*(np.logspace(np.log10(scale_ini), np.log10(scale_end), num=nscales))**2
            if kwargs.get('scale_range_is') is 'nm':
                scale_range = (3**-2)*(1./pixel_size*np.logspace(np.log10(scale_ini), np.log10(scale_end),
                                                                          num=nscales))**2
        elif scale_spacing is 'lin':
            print '\Warning in building scale range: with log spacing and arbitrary nscales, maximum in Laplacian ' \
                  'direction can give many false positive - noisy! Better use odd- spacing with a resolution.'
            if kwargs.get('scale_range_is') is 'pixel':
                scale_range = (3**-2)*(np.linspace(scale_ini, scale_end, num=nscales))**2
            if kwargs.get('scale_range_is') is 'nm':
                scale_range = (3**-2)*(1./pixel_size*np.linspace(scale_ini, scale_end, num=nscales))**2
        print '\tAnalyzing', nscales, 'scales (from t = %.2f to t = %.2f):' % (scale_range[0],
                                                                               scale_range[len(scale_range)-1])
        laplacian_mod = np.zeros(shape=(nscales, image.shape[0], image.shape[1]), dtype=float)
        l = np.zeros(shape=(nscales, image.shape[0], image.shape[1]), dtype=float)
        lx = np.zeros(shape=(nscales, image.shape[0], image.shape[1]), dtype=float)
        ly = np.zeros(shape=(nscales, image.shape[0], image.shape[1]), dtype=float)
        for n in range(1, nscales+1):
            t = scale_range[n-1]
            print '\t\tt = %.2f \t(blob diameter = %.1f pixels(analysis) = %.1f nm(physical unit))' \
                  % (t, 3*np.sqrt(t), 3*np.sqrt(t)*pixel_size)
            scalespace = compute_space_derivatives(image, t)
            laplacian = scalespace.get('lxx') + scalespace.get('lyy')
            # C = 0.5
            # laplacian_mod[n-1] = t*(lx**2 + ly**2) + C*t**2*(lxx**2 + lyy**2 + 2*lxy)  # quasi quadrature term Li98b
            laplacian_mod[n-1] = -(t ** gamma * laplacian)

            # # for descriptor, we use lin-transform Voronoi density map:
            if image_descr is not None:
                scalespace = compute_space_derivatives(image_descr, t)
                l[n-1] = scalespace.get('l'); lx[n-1] = -1. * scalespace.get('lx'); ly[n-1] = -1. * scalespace.get('ly')
        dila = np.ones(shape=(max_filter_depth, max_filter_width, max_filter_width), dtype=float)
        dila[max_filter_depth/2, max_filter_width/2, max_filter_width/2] = 0
        laplacian_mod_dila = maximum_filter(laplacian_mod, footprint=dila, mode='constant', cval=0)  # -float('inf'))
        local_maxima_locations = np.where(laplacian_mod >= laplacian_mod_dila)

        argmaxgrad = local_maxima_locations[1:]  # tuple, true values
        strength = laplacian_mod[local_maxima_locations]
        scale = scale_range[local_maxima_locations[0]]  # (1 + local_maxima_locations[0])  # 1d array

        # from src import statistics as stat
        # stat.plot_hist(strength, num_bins=50, xlabel=r'strength (Laplacian) NO THRESHOLD')
        # # stat.plot_boxplot(strength, bptype='violin', ylabel=r'strength (Laplacian) NO THRESHOLD')
        # plt.savefig('strength_boxplot.pdf', bbox_inches='tight')
        diameter = pixel_size * 3 * np.sqrt(scale)

        # begin: thresholding before orientation computation (used to be in find_feature)
        thresholding = kwargs.get('thresholding', False)
        if thresholding:
            print '\tThresholding...',
            threshold_percent = kwargs.get('threshold_percent', 0.4)  # default 40% down the maximum response
            threshold = np.max(strength) - threshold_percent * (np.max(strength) - np.min(strength))
            temp = np.where(strength > threshold)  # tuple
            argmaxgrad_threshold = (argmaxgrad[0][temp[0]], argmaxgrad[1][temp[0]])
            argmaxgrad = argmaxgrad_threshold
            scale_threshold = scale[temp[0]]
            scale = scale_threshold
            strength = strength[temp[0]]
            diameter = pixel_size * 3 * np.sqrt(scale)
            print 'Done.'
        else:
            threshold = -float('inf')
        num_features = kwargs.get('num_features', 'all')
        if num_features is not 'all':
            num_all_features = strength.shape[0]
            if num_features > num_all_features:
                print 'Warning in get_blob: size of subset of features is larger than total number of features.'
                num_features = num_all_features
            subset = np.argsort(-strength)[0:num_features]  # strongest num_features in terms of strength-laplacian
            # subset = np.random.randint(0, high=num_all_features, size=num_features, dtype='l')
            scale = scale[subset]; strength = strength[subset]; diameter = diameter[subset]
            argmaxgrad = argmaxgrad[0][subset], argmaxgrad[1][subset]
        # end

        print '\t\tNumber of detected features ', scale.shape
        aux = 0
        if compute_orientation:  # descriptors
            print '\tComputing feature orientations...'
            n_bins_ori = kwargs.get('n_bins_ori', 36)
            peak_ratio = kwargs.get('peak_ratio', 0.8)
            sigma_ori_times = kwargs.get('sigma_ori_times', 1.5)
            smooth_cycles = kwargs.get('smooth_cycles', 1)
            window_ori_radtimes = kwargs.get('window_ori_radtimes', 3)
            sigma_descr_times = kwargs.get('sigma_descr_times', 1.5)
            window_descr_radtimes = kwargs.get('window_descr_radtimes', 3)
            n_hist = kwargs.get('n_hist', 16)
            n_bins_descr = kwargs.get('n_bins_descr', 8)
            threshold_sat = kwargs.get('threshold_sat', 0.2)
            for ii, by in enumerate(argmaxgrad[1]):
                scale_index = np.where(scale_range == scale[ii])
                bx = argmaxgrad[0][ii]
                sigma_ori = sigma_ori_times * np.sqrt(scale[ii])
                # 1.5 is gaussian sigma for orientation assignment (lambda_ori in Re13)
                radius_ori = int(np.floor(window_ori_radtimes * sigma_ori))
                # radius of the (squared) region in ori assignment

                # print '\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& NEW POINT &&&&&&&&&&&'
                # print '(bx,by)=', bx, by, '\t orientation_hist \pm ', radius_ori
                # if scale[ii] > 32: plot_graphics = 1
                # else: plot_graphics = 0
                ori = orientation_hist(l[scale_index[0]], bx, by, lx[scale_index[0]][0], ly[scale_index[0]][0],
                                          radius_ori, sigma_ori, n_bins_ori, peak_ratio, smooth_cycles,
                                          image_original=image,
                                          plot_graphics=plot_graphics)
                orientation.append(ori)
                # print 'Orientation(s) computed.'
                if compute_sift_descr:
                    if ii == 0: print '\t Computing sift descriptor(s)...',
                    sigma_descr = sigma_descr_times*np.sqrt(scale[ii])
                    radius_descr = int(np.floor(window_descr_radtimes * sigma_descr))
                    # print '\n(bx,by)=', bx, by, '\t sift_descriptor \pm ', radius_descr
                    hist = sift_descriptor(l[scale_index[0]][0], bx, by, lx[scale_index[0]][0],
                                           ly[scale_index[0]][0], radius_descr, ori, n_hist, n_bins_descr,
                                           threshold_sat, plot_graphics, image_original=image)
                    histogram_descr.append(hist)
                    # if hist != []: aux = aux + len(hist)
                    if hist == []:
                        # print 'sift descriptor in Nan region -> ori set to [], even if in non-nan region.'
                        orientation[-1] = []
                    # print 'Done.'
            # print 'number of no-empty elements in histogram_descr = ', aux
            print '\t\tDONE'

    return {'argmaxgrad': argmaxgrad, 'strength': strength, 'scale': scale, 'scale_range': scale_range,
            'orientation': orientation, 'histogram_descr': histogram_descr, 'diameter': diameter}


def get_edge(image, kwargs):
    """
    From FindCluster, in fact: The factor to further reduce the analysis_pixel_size once you have a region to start
    finding clusters in e.g., if factor = 5 and analysis_pixel_size = 10nm then the program looks for clusters
    with each pixel size being 2nm params.precision: double The precision of a localization, in nm.
    This value is used as the sigma value in a 2D Gaussian point - spread function

    If nscales>1, then automatic scale selection procedure is performed.

    Input:
    ---------
    image (2d-array)
    t (float): the variance of the Gaussian kernels. Default = 1
    automatic_nscale (integer): number of scales for automatic scale selection. Loop t=2*n

    Output:
    ---------
    dictionary with
     'argmaxgrad' (tupla for 2d-array - image): argmaxgrad
     'strength' (2d-array): edgestrength
     'scale' (1d-array length argmaxgrad): selected scale parameter (variance of Gaussian kernel)
    """

    gamma = 0.5  # parameter to normalize derivatives (automatic scale selection, strongest response)
    t = kwargs.get('t')
    nscales = kwargs.get('nscales', 1)

    if nscales == 1:
        scalespace = compute_space_derivatives(image, t)
        lx = scalespace.get('lx'); ly = scalespace.get('ly')
        lxx = scalespace.get('lxx'); lxy = scalespace.get('lxy'); lyy = scalespace.get('lyy')
        lxxx = scalespace.get('lxxx'); lxxy = scalespace.get('lxxy')
        lxyy = scalespace.get('lxyy'); lyyy = scalespace.get('lyyy')
        lv = np.sqrt(lx ** 2 + ly ** 2); lvv = lx ** 2 * lxx + 2 * lx * ly * lxy + ly ** 2 * lyy
        lvvv = lx ** 3 * lxxx + 3 * lx ** 2 * ly * lxxy + 3 * lx * ly ** 2 * lxyy + ly ** 3 * lyyy
        strength = t ** (gamma * 0.5) * lv

        # find zero crossing of lvv
        lvvp = np.zeros(lvv.shape)
        lvvp[np.where(lvv > 0.)] = 1.  # greater-than-zero elements
        lvv0 = np.logical_or(abs(lvvp - np.roll(lvvp, 1, axis=1)), abs(lvvp - np.roll(lvvp, 1, axis=0)))
        lvvvn = np.zeros(lvv.shape)
        lvvvn[np.where(lvvv <= 0.)] = 1.  # >0 elements
        maxgrad = np.logical_and(lvvvn, lvv0)
        argmaxgrad = np.where(maxgrad)  # tuple, true values (=1)

        scale = np.ones(len(argmaxgrad[0])) * t

    thresholding = kwargs.get('thresholding', False)
    if thresholding:
        print '\tThresholding...',
        threshold_percent = kwargs.get('threshold_percent', 0.4)
        threshold = np.max(strength) - threshold_percent * (np.max(strength) -
                                                                   np.min(strength))

        temp = np.where(strength > threshold)  # tuple
        argmaxgrad_threshold = (argmaxgrad[0][temp[0]], argmaxgrad[1][temp[0]])
        argmaxgrad = argmaxgrad_threshold
        scale_threshold = scale[temp[0]]
        scale = scale_threshold
        strength = strength[temp[0]]
        print 'Done.'
    else:
        threshold = -float('inf')

    image_roi = np.zeros(image.shape)
    image_roi[argmaxgrad] = 1
    image_roi = flood_fill(image_roi)

    # return {'argmaxgrad': argmaxgrad, 'strength': strength, 'scale': scale}  # tuple, 2d-array and 1d-array
    return {'argmaxgrad': argmaxgrad, 'strength': strength, 'scale': scale,
            'image_roi': image_roi}


def compute_space_derivatives(image, t=10):
    """
    Computes the scale-space respresentation of an image and its derivatives.

    Input:
    --------
    image (2d- array)
    t (float): is the variance of the Gaussian kernel. Default t=10

    Output:
    --------
    dictionary {'l': l, 'lx': lx, 'ly': ly, 'lxx': lxx, 'lxy': lxy, 'lyy': lyy, 'lxxx': lxxx,
            'lxyy': lxyy, 'lxxy': lxxy, 'lyyy': lyyy}
            containing ndarray 1D with the scale-space representation and its first 3 derivatives.

    """

    from scipy import ndimage

    gauss_kernels = get_gauss_filter(t)
    g = gauss_kernels.get('g')
    dg = gauss_kernels.get('dg')
    ddg = gauss_kernels.get('ddg')
    dddg = gauss_kernels.get('dddg')

    l = ndimage.convolve(image, g[None, ...], mode='constant', cval=0.0)
    l = ndimage.convolve(l, g[..., None], mode='constant', cval=0.0)
    lx = ndimage.convolve(ndimage.filters.convolve(image, dg[None, ...], mode='constant', cval=0.0),
                          g[..., None], mode='constant', cval=0.0)
    ly = ndimage.convolve(ndimage.convolve(image, dg[..., None], mode='constant', cval=0.0),
                          g[None, ...], mode='constant', cval=0.0)
    lxx = ndimage.convolve(ndimage.convolve(image, ddg[None, ...], mode='constant', cval=0.0),
                           g[..., None], mode='constant', cval=0.0)
    lyy = ndimage.convolve(ndimage.convolve(image, ddg[..., None], mode='constant', cval=0.0),
                           g[None, ...], mode='constant', cval=0.0)
    lxy = ndimage.convolve(ndimage.convolve(image, dg[None, ...], mode='constant', cval=0.0),
                           dg[..., None], mode='constant', cval=0.0)
    lxxx = ndimage.convolve(ndimage.convolve(image, dddg[None, ...], mode='constant', cval=0.0),
                            g[..., None], mode='constant', cval=0.0)
    lyyy = ndimage.convolve(ndimage.convolve(image, dddg[..., None], mode='constant', cval=0.0),
                            g[None, ...], mode='constant', cval=0.0)
    lxyy = ndimage.convolve(ndimage.convolve(image, ddg[..., None], mode='constant', cval=0.0),
                            dg[None, ...], mode='constant', cval=0.0)
    lxxy = ndimage.convolve(ndimage.convolve(image, dg[..., None], mode='constant', cval=0.0),
                            ddg[None, ...], mode='constant', cval=0.0)

    return {'l': l, 'lx': lx, 'ly': ly, 'lxx': lxx, 'lxy': lxy, 'lyy': lyy, 'lxxx': lxxx, 'lxyy': lxyy, 'lxxy': lxxy,
            'lyyy': lyyy}


def get_gauss_filter(t=10):
    """
    Computes the Gaussian kernel and its derivatives. Used for smoothing (scale-space approach of features detections)

    Input:
    --------
    t (float): is the variance of the Gaussian kernel

    Output:
    --------
    dictionary {'g': g, 'dg': dg, 'ddg': ddg, 'dddg': dddg}
        containing ndarray 1D with the filter and the first 3 derivatives.
    """
    # x = np.linspace(-5 * np.sqrt(t), 5 * np.sqrt(t), num=np.floor(5 * np.sqrt(t) + 5 * np.sqrt(t))+1)
    t0 = t
    num_steps = 2. * round((2 * 3 * np.sqrt(t0) + 1) / 2) - 1  # with 3, 99.7% of data
    x = np.linspace(-3 * np.sqrt(t0), 3 * np.sqrt(t0), num=num_steps)  # sensitive to this, of course!
    s = np.sqrt(t)

    g = 1. / (s * np.sqrt(2 * np.pi)) * np.exp(-1 * (x * x) / (2 * s * s))
    dg = -1 * x / np.sqrt(2 * np.pi * t ** 3) * np.exp(-1 * x ** 2 / (2 * t))
    ddg = ((x - np.sqrt(t)) * (x + np.sqrt(t))) / (np.sqrt(2 * np.pi * t ** 5)) * np.exp(-x ** 2 / (2 * t))
    dddg = ((-x * (x ** 2 - 3 * t)) / np.sqrt(2 * np.pi * t ** 7)) * np.exp(-x ** 2 / (2 * t))

    # plt.plot(x, g, 'k:', label='$t=\sigma^2=20$')
    # plt.xlabel('pixel')
    # plt.legend(loc='upper left')
    # plt.xlim(-20, 20); plt.ylim(0, 0.3)
    # plt.show()
    # plt.hold(True)

    return {'g': g, 'dg': dg, 'ddg': ddg, 'dddg': dddg}


def plot_feature(image, feature, feature_name='blob', cmap='gray', interpolation='none', norm=None,
                 plot_axis='on', blob_color='strength', ori_color=None, ori_cmap=None):
    """
    Function for plotting scale-space feature

    Input:
    --------
    image (ndarray)
    feature (dictionary)
    cmap: map of the image
    norm (string): 'log', 'lin'. Scale of imshow. E.g.,
        if norm is 'log':  plt.imshow(image.T, interpolation='none', cmap='jet', origin='lower', norm=LogNorm())
    blob_color = 'strength', 'scale', 'class'
    ori_color (list) = if not None, for all orientations, assigned classes
    ori_cmap = if not None, discrete cmap for vocabulary - output clustering/unsupervised learing.

    Output:
    --------

    """
    from matplotlib.colors import LogNorm
    import matplotlib.colors as colors
    from src import utilities as util
    from collections import Counter

    fig, ax = plot_image(image, cmap=cmap, interpolation=interpolation, norm=norm, plot_axis=plot_axis, hold=True)
    plt.hold(1)
    argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    scale = feature.get('scale')  # 1d-array
    strength = feature.get('strength')
    orientation = feature.get('orientation', [])

    if ori_color is not None and ori_cmap is None: ori_cmap = plt.cm.get_cmap('Set1')
    if blob_color == 'strength' or blob_color == 'scale':
        values = []
        if blob_color is 'strength': values = strength
        elif blob_color is 'scale': values = np.sqrt(scale) * 2*1.5
        cnorm = colors.Normalize(vmin=np.min(values), vmax=1.05*np.max(values))
        blob_cmap = plt.cm.ScalarMappable(norm=cnorm, cmap=plt.cm.gray)
        blob_cmap.set_array(values)
    elif blob_color is 'class':
        if ori_cmap is None: raise ValueError('Error in iproc.plot_feature: introduce ori_color.')
        blob_cmap = ori_cmap
    else: ValueError('Error in iproc.plot: introduce blob_color \in (strength, scale, class)')

    if feature_name == 'edge':
        plt.plot(argmaxgrad[0], argmaxgrad[1], 'r.', markersize=5)

    if feature_name == 'blob':
        cval = (2 * np.pi * np.linspace(0, 1, num=50))
        ucirc = np.array([np.cos(cval), np.sin(cval)])
        hist_ind = 0   # plot orientation with colorcode from clustering algorithm (label)
        for ii, by in enumerate(argmaxgrad[1]):
            bx = argmaxgrad[0][ii]
            # # plot arrows dominant orientations
            if len(orientation) > 0:  # if descriptors have been computed
                mean = []
                # if orientation[ii] == []: print 'empty at hist_ind = ', hist_ind, ' ; argmaxgrad ii = ', ii
                for jj, ori in enumerate(orientation[ii]):  # if [] => skip, this blob is not in the analysis
                    if ori_color is not None:  # colors according to unsupervised
                        o_color = ori_cmap(ori_color[hist_ind])
                        mean.append(ori_color[hist_ind])
                        hist_ind += 1
                    elif blob_color == 'strength': o_color = blob_cmap.to_rgba(strength[ii])
                    elif blob_color == 'scale': o_color = blob_cmap.to_rgba(np.sqrt(scale[ii]) * 2 * 1.5)
                    plt.plot([bx, bx + np.sqrt(scale[ii]) * 1.5 * np.cos(ori)],
                             [by, by + np.sqrt(scale[ii]) * 1.5 * np.sin(ori)],
                              color=o_color, linewidth=1)
                if len(orientation[ii]) > 0 and ori_color is not None: mean = Counter(mean).most_common(1)[0][0]
            # # plot blobs - detected features
            if blob_color == 'strength': b_color = blob_cmap.to_rgba(strength[ii])
            elif blob_color == 'scale': b_color = blob_cmap.to_rgba(np.sqrt(scale[ii]) * 2 * 1.5)
            elif blob_color == 'class':
                if len(orientation[ii]) == 0: b_color = 'None'
                else: b_color = blob_cmap(mean)
            ax.plot(ucirc[0, :]*np.sqrt(scale[ii])*1.5 + bx, ucirc[1, :]*np.sqrt(scale[ii])*1.5 +
                    by, color=b_color, linewidth=1)
            # ax.text(bx, by, '(%d,%d; %.1f; %.2f' % (bx, by, 3 * np.sqrt(scale[ii]), strength), color='r')
        # fig.colorbar(blob_cmap, label='max$_t\{\Delta_{\gamma-norm}\}$')
        # fig.colorbar(blob_cmap, label='3$\sqrt{t}$ [pixel - analysis]')

    ax.set_xlim(0, image.T.shape[1]); ax.set_ylim(0, image.T.shape[0])
    ax.set_aspect('equal', adjustable='box')

    return ax.figure


def orientation_hist(image, bx, by, lx, ly, radius, sigma_ori, n_bins_ori, peak_ratio=0.8, smooth_cycles=2,
                     image_original=None, plot_graphics=True):
    """
    Computes a gradient orientation histogram at a specified pixel
    compute gradient values, orientations and the weights over the pixel neighborhood
    """

    from matplotlib.colors import LogNorm

    hist = np.zeros(shape=n_bins_ori)  # lx.shape=image.shape (row,cols), transpose from reading image =>
    # bx-by=row-col in this case!
    grad_x = lx[bx-radius:bx+radius+1, by-radius:by+radius+1]
    grad_y = ly[bx-radius:bx+radius+1, by-radius:by+radius+1]
    # print 'lx.shape, ly.shape =', lx.shape, ly.shape
    # fig, ax = plot_image(lx.T, cmap='jet', interpolation='none', hold=True)
    if np.isnan(grad_x).any() or np.isnan(grad_y).any() or \
       bx-radius < 0 or bx+radius > image.shape[1] or by-radius < 0 or by + radius > image.shape[2]:
        print 'out'
        return []
    locgrid_x, locgrid_y = np.meshgrid(np.linspace(bx - radius, bx + radius, 2*radius+1), # cart. ref. frame
                                       np.linspace(by - radius, by + radius, 2*radius + 1))
    weights = np.exp(-((locgrid_x-bx)**2 + (locgrid_y-by)**2)/(2*sigma_ori**2))
    ori = np.arctan2(grad_x, grad_y)  # in the image.T reference (rotate) -> (bx,by)
    contr = np.multiply(weights, np.sqrt(grad_x**2 + grad_y**2))  # contribution to the histogram

    # build histogram
    for k in range((2*radius+1)*(2*radius+1)):
        i = int(np.round(ori[k/int(2*radius+1)][k % int(2*radius+1)]/(2*np.pi/n_bins_ori)))
        hist[i] += contr[k/int(2*radius+1)][k % int(2*radius+1)]  # in python hist[-1] = hist[last position], etc.
    # smooth histogram
    hist = smooth_hist(hist, smooth_cycles)
    # find maxima
    argpeak = np.argmax(hist)
    from scipy.ndimage.filters import maximum_filter
    max_filter_width = 5; dila = np.ones(shape=max_filter_width, dtype=float); dila[max_filter_width/2] = 0
    argmax_all = np.where(hist > maximum_filter(hist, footprint=dila, mode='wrap'))
    argmax = argmax_all[0][hist[argmax_all] >= peak_ratio*hist[argpeak]]

    if plot_graphics:
        # plots
        plt.figure()
        ax = plt.subplot(1, 2, 1)
        plt.imshow(image_original.T, interpolation='none', cmap='gray', origin='lower')  #, norm=LogNorm())
        plt.hold(True)
        for i in range(0, 2*radius+1, 2):
            y = by + i - radius
            for j in range(0, 2*radius+1, 2):
                x = bx + j - radius  # col
                ax.arrow(x, y, 10000*weights[j, i]*grad_y[j, i], 10000*weights[j, i]*grad_x[j, i],
                         head_width=0.2, head_length=0.2, fc='r', ec='r', width=0.1)
        plt.show()
        plt.subplot(1, 2, 2); plt.bar(np.arange(n_bins_ori), hist, color='0.5', align='center')
        plt.xticks(np.linspace(0, n_bins_ori, 5), ['0', '\pi/2',  '\pi',  '3\pi/2', '2\pi'])
        plt.axhline(peak_ratio*hist[argpeak], color='k', ls=':')
        plt.show(); plt.hold(True)
        plt.plot(argmax, hist[argmax], 'r*', markersize=10)

    return argmax * 2 * np.pi / n_bins_ori


def smooth_hist(hist, smooth_cycles=6):
    """
    Function that smoothes the histogram according to (Rey-Otero et al., 2014, pg. 386).
    That is: a number of times a circular convolution with filter [1 1 1].1/3
    """
    bins = hist.shape[0]
    for ii in range(smooth_cycles*hist.shape[0]):
        hist[ii % bins] = 1./3*(hist[(ii-1) % bins]+hist[ii % bins]+hist[(ii+1) % bins])

    return hist


def sift_descriptor(image, bx, by, lx, ly, radius, orientation, n_hist=16, n_bins_descr=8,
                    threshold_sat=0.2, plot_graphics=True, image_original=None):
    """
    Calculate sample's histogram array coords rotated relative to ori.

    Output:
    -----------
    """
    from matplotlib.colors import LogNorm

    orientation = np.asarray(orientation)

    if orientation.size == 0:
        # print '=============================== Because no main orienations, compute_orientation returns []'
        return []

    # image = image_test
    # plt.close('all')
    # n_hist = 16
    # n_bins_descr = 8
    # window_descr_radtimes = 4
    # argmaxgrad = feature.get('argmaxgrad')
    # scale = feature.get('scale')  # 1d-array
    # scale_range = feature.get('scale_range')
    # orientation = feature.get('orientation')
    # sigma_descr_times = 1.5
    # sigma_descr = sigma_descr_times * np.sqrt(scale[0])
    # radius = int(np.floor(window_descr_radtimes * sigma_descr))
    # bx = argmaxgrad[0][0]
    # by = argmaxgrad[1][0]
    # ori = orientation[0][0]  #3*np.pi/4.
    #
    # scalespace = iproc.compute_space_derivatives(image, sigma)
    # l = scalespace.get('l')
    # lx= -1.*scalespace.get('lx')  # local minima, gradient descent, pointing from white to black (1->0)
    # ly = -1.*scalespace.get('ly')

    n_hist_x = np.sqrt(n_hist)
    grid_x, grid_y = np.meshgrid(np.arange(0, image.shape[0]), np.arange(0, image.shape[1]))  # OK, 'image' is
    # rotated: #rows=xlim; #cols=ylim

    histogram = []
    for ori in orientation:
        R = np.array([[np.cos(ori), -np.sin(ori)],
                      [np.sin(ori), np.cos(ori)]])
        R_inv = R.T
        corners = (R.dot(radius*np.array([[-1, -1, 1, 1], [-1, 1, 1, -1]]))+np.array([[bx], [by]])).T
        m = (corners[[1, 2, 3, 0], 1] - corners[..., 1])/(corners[[1, 2, 3, 0], 0] - corners[..., 0])
        c = corners[..., 1]-corners[..., 0]*m  # notes 18/5/17
        ind = int((ori-np.pi/2.)//(0.5*np.pi))
        ind0 = ind % 4; ind1 = (ind+1) % 4; ind2 = (ind+2) % 4; ind3 = (ind+3) % 4
        if ori/(np.pi*0.5) % 1 == 0:
            pos = np.where((grid_x < corners[ind0, 0])
                           & (grid_y > m[(ind1-1) % 4]*grid_x + c[(ind1-1) % 4])
                           & (grid_x > corners[(ind2-1) % 4, 0])
                           & (grid_y < m[(ind3-1) % 4]*grid_x + c[(ind3-1) % 4]))
        else:
            pos = np.where((grid_y > m[ind0] * grid_x + c[ind0])
                           & (grid_y > m[ind1] * grid_x + c[ind1])
                           & (grid_y < m[ind2] * grid_x + c[ind2])
                           & (grid_y < m[ind3] * grid_x + c[ind3]))

        pos_ref = R_inv.dot(np.array([pos[1]-bx, pos[0]-by]))
        hist_ind = np.floor((pos_ref[0]+radius)//(2.*radius/n_hist_x*1.000001) + n_hist_x*((pos_ref[1]+radius)//(
                   2.*radius/n_hist_x)*1.000001))
        grad_x = lx[pos[1], pos[0]]; grad_y = ly[pos[1], pos[0]]
        if np.isnan(grad_x).any() or np.isnan(grad_y).any():
            # print '=============================== blob in nan region. compute_orientation returns []'
            return []
        ori_all = (np.arctan2(grad_x, grad_y) - ori)*180.//np.pi % 360

        weights = np.exp(-((pos[0]-by)**2 + (pos[1]-bx)**2)/(2*radius**2))
        contr = np.multiply(weights, np.sqrt(grad_x**2 + grad_y**2))  # contribution to the histogram

        # build the histogram
        hist = np.zeros(shape=(n_hist*n_bins_descr))
        ori_all_bins = np.floor(ori_all / (360 / n_bins_descr))
        # print '\n(bx,by)=', bx, by, ' ******************ori=', ori
        # print np.min(hist_ind), np.min(np.arctan2(grad_x, grad_y) - ori), np.min(ori_all), np.min(ori_all_bins)
        # print np.max(hist_ind), np.max(np.arctan2(grad_x, grad_y) - ori), np.max(ori_all), np.max(ori_all_bins)
        # print 'hist_ind[0]=', hist_ind[0]
        # print ori_all[0], ori_all_bins[0]
        bins_temp = (n_bins_descr * hist_ind + ori_all_bins).astype(int)
        for i, bins in enumerate(bins_temp):
            hist[bins] += contr[i]

        # normalize descriptor vectors
        hist /= np.linalg.norm(hist)
        hist[np.where(hist > threshold_sat)] = threshold_sat  # we reduce the influence of large gradient
        # magnitudes by thresholding the values in the unit feature vector to each be no larger than 0.2,
        # and then renormalizing to unit length.
        hist /= np.linalg.norm(hist)

        # print bx, by
        histogram.append(hist)

        ind = 4 * (3 - 0 // 4) + 0 % 4
        if plot_graphics:
            if image_original is not None:
                image = image_original

            # plt.figure();  plt.plot(grid_x, grid_y, 'r.', linestyle='none')
            # plt.plot(corners[..., 0], corners[..., 1], 'k.', linestyle='none')
            # plt.plot(grid_x[pos], grid_y[pos], 'k.', linestyle='none', alpha=1); plt.axis('equal')
            # plt.text(corners[0, 0]+2, corners[0, 1]+2, 'zero')
            # plt.text(corners[1, 0]+2, corners[1, 1]+2, 'one')
            # plt.text(corners[2, 0]+2, corners[2, 1]+2, 'two')
            # plt.text(corners[3, 0]+2, corners[3, 1]+2, 'three')
            # plt.axis('equal')

            # plt.figure(); plt.imshow(image.T, interpolation='none', cmap='gray', origin='low'); ax = plt.axes()
            # ax.autoscale(False); ax.hold(True)
            # sc = plt.scatter(pos[1], pos[0], c=hist_ind, alpha=0.4); plt.hold(True)
            # plt.colorbar(sc)
            # plt.text(corners[0, 0]+2, corners[0, 1]+2, '0', color='w')
            # plt.text(corners[1, 0]+2, corners[1, 1]+2, '1', color='w')
            # plt.text(corners[2, 0]+2, corners[2, 1]+2, '2', color='w')
            # plt.text(corners[3, 0]+2, corners[3, 1]+2, '3', color='w')
            # plt.show()

            plt.figure(); plt.title('center = (%d,%d), ori = %f' % (bx, by, ori))
            ax = plt.subplot(1, 2, 1)
            plt.imshow(image.T, interpolation='none', cmap='gray', origin='lower')  # , norm=LogNorm())
            # plt.imshow(image.T, interpolation='none', cmap='gray', origin='lower')
            plt.show(); plt.hold(True)
            ax.autoscale(False); ax.hold(True)

            sc = plt.scatter(pos[1], pos[0], c=ori_all_bins, alpha=0.4)
            cbar = plt.colorbar(sc, ticks=np.linspace(0, n_bins_descr, 5), aspect=50)
            cbar.ax.set_yticklabels(['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
            # cbar.ax.set_ylabel('relative orientation [rad]'); plt.show()
            step = 2./n_hist_x
            for ii in np.arange(0, n_hist_x+1):
                xy0 = (R.dot(radius*np.array([-1+step*ii, -1]).T+0*np.array([0, 1]).T)+np.array([bx, by]).T).T
                xy1 = (R.dot(radius*np.array([-1+step*ii, -1]).T+2*radius*np.array([0, 1]).T)+np.array([bx, by]).T).T
                plt.plot([xy0[0], xy1[0]], [xy0[1], xy1[1]], 'r-')
                plt.hold(True)
                xy0 = (R.dot(radius*np.array([-1, -1+step*ii]).T+0*np.array([1, 0]).T)+np.array([bx, by]).T).T
                xy1 = (R.dot(radius*np.array([-1, -1+step*ii]).T+2*radius*np.array([1, 0]).T)+np.array([bx, by]).T).T
                plt.plot([xy0[0], xy1[0]], [xy0[1], xy1[1]], 'r-')
            ax.arrow(bx, by, 15*np.cos(ori), 15*np.sin(ori), head_width=3, head_length=3, fc='w', ec='w', alpha=1,
                     linewidth=2)
            plt.show()

            max_ori = np.max(hist)
            for ii in range(0, n_hist):
                ind = 4*(3-ii//4) + ii % 4
                plt.subplot(4, 8, (ii//4+1)*4+ii+1)
                # plt.subplot(4, 4, ii+1)
                # print 'checking histogram sift:'
                # print n_bins_descr*ind, ' to ', n_bins_descr*ind+n_bins_descr
                plt.bar(np.arange(n_bins_descr), hist[n_bins_descr*ind:n_bins_descr*ind+n_bins_descr], color='0.5',
                        align='center')
                plt.xticks(np.linspace(0, n_bins_descr, 5), ['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
                plt.ylim([0, max_ori])
                plt.hold(True)

    return histogram

