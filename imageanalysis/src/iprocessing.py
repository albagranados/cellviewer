import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# configure latex plots
matplotlib.rcParams['text.usetex'] = True  # uselatex labels
matplotlib.rcParams['font.family'] = 'serif'


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

    area_region = range_x*range_y  # input total image area
    if pixel_size is None:
        scale_pixel_size = math.sqrt(2) * area_region / float(num_points)
    # average_pointdensity = float(num_points)/area_region

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

    plt.show()
    if plot_axis is 'off':
        plt.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    plt.hold(hold)
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


def find_feature(image, kwargs):
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

    Output:
    ---------
    edgemat (dictionary)

    """
    t = kwargs.get('t', 10)
    feature_name = kwargs.get('feature_name', 'blob')
    thresholding = kwargs.get('thresholding', False)
    automatic_nscales = kwargs.get('automatic_nscales', 1)

    if feature_name == 'edge':
        feature = get_edge(image, t, automatic_nscales)

    if feature_name == 'blob':
        feature = get_blob(image, kwargs)

    #  # begin: this code below has been moved to get_blob -> we compute orientations just in the thresholded features.
    # argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    # featurestrength = feature.get('featurestrength')  # 2d-array
    # tnew = feature.get('tnew')
    # if thresholding:
    #     threshold_percent = kwargs.get('threshold_percent', 0.4)  # default 40% down the maximum response
    #     threshold = np.max(featurestrength[argmaxgrad])-threshold_percent*(np.max(featurestrength[argmaxgrad]) -
    #                                                                        np.min(featurestrength[argmaxgrad]))
    # else:
    #     threshold = -float('inf')
    #
    # feature['threshold'] = threshold
    #
    # temp = np.where(featurestrength[argmaxgrad] > threshold)  # tuple
    # argmaxgrad_threshold = (argmaxgrad[0][temp[0]], argmaxgrad[1][temp[0]])
    # feature['argmaxgrad'] = argmaxgrad_threshold
    # tnew_threshold = tnew[temp[0]]
    # feature['tnew'] = tnew_threshold
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

    print('Done.')
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


def get_blob(image, kwargs):
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
     'featurestrength' (2d-array): blobstrength
     'tnew' (1d-array length argmaxgrad): selected scale parameter (variance of Gaussian kernel)

    """
    from scipy.ndimage.filters import maximum_filter

    t = kwargs.get('t', 1)
    original_pixel_size = kwargs.get('original_pixel_size', 1)
    scale_pixel_size = kwargs.get('scale_pixel_size', 1)
    pixel_size = scale_pixel_size * original_pixel_size
    resolution = kwargs.get('resolution', 1)
    automatic_nscales = kwargs.get('automatic_nscales', 1)
    max_filter_width = kwargs.get('max_filter_width', 3)
    max_filter_depth = kwargs.get('max_filter_depth', max_filter_width)
    scale_spacing = kwargs.get('scale_spacing', 'log')
    scale_ini = kwargs.get('scale_ini', 1)
    scale_end = kwargs.get('scale_end', 10)
    compute_orientation = kwargs.get('compute_orientation', False); orientation = []
    compute_sift_descr = kwargs.get('compute_sift_descr', False); histogram_descr = []
    plot_graphics = kwargs.get('plot_graphics', False)

    gamma = 1  # parameter normalized derivatives
    if automatic_nscales == 1:
        scalespace = compute_space_derivatives(image, t)
        lxx = scalespace.get('lxx')
        lyy = scalespace.get('lyy')
        laplacian = lxx + lyy
        featurestrength = -t ** gamma * laplacian  # +t** -> find minima

        dila = np.ones(shape=(max_filter_width, max_filter_width), dtype=float)
        dila[max_filter_width/2, max_filter_width/2] = 0
        featurestrength_dila = maximum_filter(featurestrength, footprint=dila, mode='constant', cval=-float('inf'))
        argmaxgrad = np.where(featurestrength > featurestrength_dila)

        tnew = np.ones(len(argmaxgrad[0])) * t
        scale_range = tnew
    else:
        if scale_spacing is 'odd':
            if kwargs.get('scale_range_is') is 'nm':
                scale_ini = np.ceil((scale_ini / pixel_size * pixel_size / resolution - 1) / 2.)
                scale_end = np.ceil((scale_end / pixel_size * pixel_size / resolution - 1) / 2.)
                automatic_nscales = int(scale_end - scale_ini) + 1
                max_filter_depth = automatic_nscales + 1
            elif kwargs.get('scale_range_is') is 'pixel':
                scale_ini = np.ceil((scale_ini - 1) / 2.)
                scale_end = np.ceil((scale_end - 1) / 2.)
                automatic_nscales = int(scale_end - scale_ini) + 1
                max_filter_depth = automatic_nscales + 1
            scale_range = (3 ** -2) * (resolution / pixel_size * (2 * np.linspace(scale_ini, scale_end,
                                                                                  num=automatic_nscales) + 1)) ** 2
        elif scale_spacing is 'log':
            if kwargs.get('scale_range_is') is 'pixel':
                scale_range = (3**-2)*(np.logspace(np.log10(scale_ini), np.log10(scale_end), num=automatic_nscales))**2
            if kwargs.get('scale_range_is') is 'nm':
                scale_range = (3**-2)*(resolution/pixel_size*np.logspace(np.log10(scale_ini), np.log10(scale_end),
                                                                          num=automatic_nscales))**2
        elif scale_spacing is 'lin':
            if kwargs.get('scale_range_is') is 'pixel':
                scale_range = (3**-2)*(np.linspace(scale_ini, scale_end, num=automatic_nscales))**2
            if kwargs.get('scale_range_is') is 'nm':
                scale_range = (3**-2)*(resolution/pixel_size*np.linspace(scale_ini, scale_end,
                                                                     num=automatic_nscales))**2

        print '\tAnalyzing', automatic_nscales, 'scales (from t=%.2f to t=%.2f):\n' %(scale_range[0], scale_range[1])
        strength = np.zeros(shape=(automatic_nscales, image.shape[0], image.shape[1]), dtype=float)
        l = np.zeros(shape=(automatic_nscales, image.shape[0], image.shape[1]), dtype=float)
        lx = np.zeros(shape=(automatic_nscales, image.shape[0], image.shape[1]), dtype=float)
        ly = np.zeros(shape=(automatic_nscales, image.shape[0], image.shape[1]), dtype=float)
        for n in range(1, automatic_nscales+1):
            t = scale_range[n-1]
            print '\t\tt = %.2f (blob diameter = %.1f pixels(analysis) = %.1f nm(pysical units))' % (t, 3*np.sqrt(t),
                                                                            3*np.sqrt(t)*pixel_size)
            scalespace = compute_space_derivatives(image, t)
            lxx = scalespace.get('lxx')
            lyy = scalespace.get('lyy')
            l[n-1] = scalespace.get('l')
            lx[n-1] = -1. * scalespace.get('lx')  # local minima, gradient descent, pointing from white to black (1->0)
            ly[n-1] = -1. * scalespace.get('ly')
            laplacian = lxx + lyy
            # C = 0.5
            # strength[n-1] = t*(lx**2 + ly**2) + C*t**2*(lxx**2 + lyy**2 + 2*lxy)  # quasi quadrature term Li98b
            strength[n-1] = -(t ** gamma * laplacian)

        dila = np.ones(shape=(max_filter_depth, max_filter_width, max_filter_width), dtype=float)
        dila[max_filter_depth/2, max_filter_width/2, max_filter_width/2] = 0
        featurestrength_dila = maximum_filter(strength, footprint=dila, mode='constant', cval=0)  # -float('inf'))
        local_maxima_locations = np.where(strength > featurestrength_dila)

        argmaxgrad = local_maxima_locations[1:]  # tuple, true values
        featurestrength = np.zeros(shape=image.shape)  # 2d array
        featurestrength[argmaxgrad] = strength[local_maxima_locations]

        tnew = scale_range[local_maxima_locations[0]]  # (1 + local_maxima_locations[0])  # 1d array

        # begin: thresholding before orientation computation (used to be in find_feature)
        thresholding = kwargs.get('thresholding', False)
        if thresholding:
            print '\tThresholding...',
            threshold_percent = kwargs.get('threshold_percent', 0.4)  # default 40% down the maximum response
            threshold = np.max(featurestrength[argmaxgrad]) - threshold_percent * (np.max(featurestrength[argmaxgrad]) -
                                                                                   np.min(featurestrength[argmaxgrad]))
            print 'Done.'
        else:
            threshold = -float('inf')

        temp = np.where(featurestrength[argmaxgrad] > threshold)  # tuple
        argmaxgrad_threshold = (argmaxgrad[0][temp[0]], argmaxgrad[1][temp[0]])
        argmaxgrad = argmaxgrad_threshold
        tnew_threshold = tnew[temp[0]]
        tnew = tnew_threshold
        image_roi = np.zeros(image.shape)
        image_roi[argmaxgrad_threshold] = 1
        image_roi = flood_fill(image_roi)
        # end

        if compute_orientation:
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
                scale_index = np.where(scale_range == tnew[ii])
                bx = argmaxgrad[0][ii]
                print '\t\t(bx,by)=(%d,%d)... ' % (bx,by),
                sigma_ori = sigma_ori_times * np.sqrt(tnew[ii])
                # 1.5 is gaussian sigma for orientation assignment (lambda_ori in Re13)
                radius = int(np.floor(window_ori_radtimes * sigma_ori))
                # radius of the (squared) region in ori assignment

                ori = orientation_hist(l[scale_index[0]], bx, by, lx[scale_index[0]][0], ly[scale_index[0]][0],
                                       radius, sigma_ori, n_bins_ori, peak_ratio, smooth_cycles, image_original=image,
                                       plot_graphics=plot_graphics)
                orientation.append(ori)
                print 'Orientation(s) computed.'
                if compute_sift_descr:
                    print '\t\t\t Computing sift descriptor(s)...',
                    sigma_descr = sigma_descr_times*np.sqrt(tnew[ii])
                    radius_descr = int(np.floor(window_descr_radtimes * sigma_descr))
                    hist = sift_descriptor(l[scale_index[0]][0], bx, by, lx[scale_index[0]][0],
                                           ly[scale_index[0]][0], radius_descr, ori, n_hist, n_bins_descr,
                                           threshold_sat, plot_graphics, image_original=image)
                    histogram_descr.append(hist)
                    print 'Done.'
            print '\t\tDONE'

    # return {'argmaxgrad': argmaxgrad, 'featurestrength': featurestrength, 'tnew': tnew, 'scale_range': scale_range,
    #         'orientation': orientation, 'histogram_descr': histogram_descr}
    return {'argmaxgrad': argmaxgrad, 'featurestrength': featurestrength, 'tnew': tnew,
            'scale_range': scale_range,
            'orientation': orientation, 'histogram_descr': histogram_descr, 'threshold': threshold,
            'image_roi': image_roi}


def get_edge(image, t=1, automatic_nscales=1):
    """
    From FindCluster, in fact: The factor to further reduce the analysis_pixel_size once you have a region to start
    finding clusters in e.g., if factor = 5 and analysis_pixel_size = 10nm then the program looks for clusters
    with each pixel size being 2nm params.precision: double The precision of a localization, in nm.
    This value is used as the sigma value in a 2D Gaussian point - spread function

    If automatic_nscales>1, then automatic scale selection procedure is performed.

    Input:
    ---------
    image (2d-array)
    t (float): the variance of the Gaussian kernels. Default = 1
    automatic_nscale (integer): number of scales for automatic scale selection. Loop t=2*n

    Output:
    ---------
    dictionary with
     'argmaxgrad' (tupla for 2d-array - image): argmaxgrad
     'featurestrength' (2d-array): edgestrength
     'tnew' (1d-array length argmaxgrad): selected scale parameter (variance of Gaussian kernel)
    """

    gamma = 0.5  # parameter to normalize derivatives (automatic scale selection, strongest response)

    if automatic_nscales == 1:
        scalespace = compute_space_derivatives(image, t)
        lx = scalespace.get('lx'); ly = scalespace.get('ly')
        lxx = scalespace.get('lxx'); lxy = scalespace.get('lxy'); lyy = scalespace.get('lyy')
        lxxx = scalespace.get('lxxx'); lxxy = scalespace.get('lxxy')
        lxyy = scalespace.get('lxyy'); lyyy = scalespace.get('lyyy')
        lv = np.sqrt(lx ** 2 + ly ** 2); lvv = lx ** 2 * lxx + 2 * lx * ly * lxy + ly ** 2 * lyy
        lvvv = lx ** 3 * lxxx + 3 * lx ** 2 * ly * lxxy + 3 * lx * ly ** 2 * lxyy + ly ** 3 * lyyy
        featurestrength = t ** (gamma * 0.5) * lv

        # find zero crossing of lvv
        lvvp = np.zeros(lvv.shape)
        lvvp[np.where(lvv > 0.)] = 1.  # greater-than-zero elements
        lvv0 = np.logical_or(abs(lvvp - np.roll(lvvp, 1, axis=1)), abs(lvvp - np.roll(lvvp, 1, axis=0)))
        lvvvn = np.zeros(lvv.shape)
        lvvvn[np.where(lvvv <= 0.)] = 1.  # >0 elements
        maxgrad = np.logical_and(lvvvn, lvv0)
        argmaxgrad = np.where(maxgrad)  # tuple, true values (=1)

        tnew = np.ones(len(argmaxgrad[0])) * t

    return {'argmaxgrad': argmaxgrad, 'featurestrength': featurestrength, 'tnew': tnew}  # tuple, 2d-array and 1d-array


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
    num_steps = 2. * round((2 * 5 * np.sqrt(t0) + 1) / 2) - 1
    x = np.linspace(-5 * np.sqrt(t0), 5 * np.sqrt(t0), num=num_steps)
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


def plot_feature(image, feature, cmap='gray', interpolation='none', norm=None, plot_axis='on',
                 feature_name='blob', blob_color=None, cluster_color=None):
    """
    Function for plotting scale-space feature

    Input:
    --------
    image (ndarray)
    feature (dictionary): {'argmaxgrad': argmaxgrad, 'featurestrength': featurestrength, 'tnew', 'scale_range'}
    cmap: map of the image
    norm (string): 'log', 'lin'. Scale of imshow. E.g.,
        if norm is 'log':
        plt.imshow(image.T, interpolation='none', cmap='jet', origin='lower', norm=LogNorm())

    Output:
    --------

    """
    from matplotlib.colors import LogNorm
    import matplotlib.colors as colors

    fig, ax = plot_image(image, cmap=cmap, interpolation=interpolation, norm=norm, plot_axis=plot_axis, hold=True)

    argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    tnew = feature.get('tnew')  # 1d-array
    featurestrength = feature.get('featurestrength')
    orientation = feature.get('orientation', [])

    if blob_color is None:
        cmap = plt.cm.gray
        values = featurestrength[argmaxgrad[0], argmaxgrad[1]]  # range(argmaxgrad[0].shape[0])
        # values = np.sqrt(tnew) * 2*1.5
        cnorm = colors.Normalize(vmin=np.min(values), vmax=np.max(values)) #colors.Normalize(vmin=0, vmax=np.max(values))
        #
        # vmax=values[-1]) . LogNorm, Normalize
        scalarmap = plt.cm.ScalarMappable(norm=cnorm, cmap=cmap)  # norm=cNorm
        scalarmap.set_array(values)
    if cluster_color is not None:
        cmap = plt.get_cmap('Set1', np.max(cluster_color) - np.min(cluster_color) + 1)
        scalarmap_clusters = plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=np.max(cluster_color)),
                                                cmap=cmap)
        scalarmap_clusters.set_array(cluster_color)

    if feature_name == 'edge':
        plt.plot(argmaxgrad[0], argmaxgrad[1], 'r.', markersize=5)

    if feature_name == 'blob':
        cval = (2 * np.pi * np.linspace(0, 1, num=50))
        ucirc = np.array([np.cos(cval), np.sin(cval)])
        hist_ind = 0   # plot orientation with colorcode from clustering algorithm (label)
        for ii, by in enumerate(argmaxgrad[1]):
            # if tnew[ii] == feature.get('scale_range')[0]:
            # if nloc is not None and nloc[ii]==4:
                bx = argmaxgrad[0][ii]
                # # if strength < 1:
                if scalarmap is not None:
                    # strength = np.sqrt(tnew[ii]) * 2 * 1.5
                    strength = featurestrength[argmaxgrad[0][ii], argmaxgrad[1][ii]]
                    blob_color = scalarmap.to_rgba(strength)  # values[ii])
                ax.plot(ucirc[0, :] * np.sqrt(tnew[ii]) * 1*1.5 + bx, ucirc[1, :] * np.sqrt(tnew[ii]) * 1*1.5 +
                        by, color=blob_color, linewidth=1.5)
                # ax.text(1.1* np.sqrt(tnew[ii]) * 1*1.5 + bx, 1.1* np.sqrt(tnew[ii]) * 1*1.5 + by, '%.2f'%strength,
                #         color='k')
                if len(orientation) > 0:
                    for jj in orientation[ii]:
                        if cluster_color is not None:
                            ori_color = scalarmap_clusters.to_rgba(cluster_color[hist_ind])
                            hist_ind += 1
                        else:
                            ori_color = blob_color
                        plt.arrow(bx, by,
                                  np.sqrt(tnew[ii]) * 1*1.5 * np.cos(jj), np.sqrt(tnew[ii]) * 1*1.5 * np.sin(jj),
                                  head_width=0, head_length=0, fc=ori_color, ec=ori_color, fill=True, width=0.2)
        # fig.colorbar(scalarmap, label='max$_t\{\Delta_{\gamma-norm}\}$')
        # fig.colorbar(scalarmap, label='3$\sqrt{t}$')

    ax.set_xlim(0, image.T.shape[1])
    ax.set_ylim(0, image.T.shape[0])
    ax.set_aspect('equal', adjustable='box')

    return ax.figure


# def compute_sift(image, feature, n_bins_ori=18, peak_ratio=0.8, smooth_cycles=1, sigma_ori_times=2,
#                  window_ori_radtimes=3):
#     """
#     Compute Scale-Invariant Feature Transformation (SIFT, "Distinctive Image Features from Scale-Invariant
#     Keypoints", Lowe, 2004)
#     [X, Y, S, TH, D]. computes the SIFT descriptors as well. Each column of D is the descriptor of the corresponding
#     frame.A descriptor is a 128 - dimensional vector. X, Y is the center of the frame, S is the scale and TH is the
#     orientation (in radians).
#
#     Input:
#     -----------
#     n (int): num bins
#
#     Output:
#     ------------
#
#     """
#     argmaxgrad = feature.get('argmaxgrad')
#     tnew = feature.get('tnew')  # 1d-array
#     scale_range = feature.get('scale_range')
#
#     l = np.zeros(shape=(len(scale_range), image.shape[0], image.shape[1]))
#     lx = np.zeros(shape=(len(scale_range), image.shape[0], image.shape[1]))
#     ly = np.zeros(shape=(len(scale_range), image.shape[0], image.shape[1]))
#     for ii in range(len(scale_range)):
#         scalespace = compute_space_derivatives(image, scale_range[ii])
#         l[ii] = scalespace.get('l')
#         lx[ii] = -1.*scalespace.get('lx')  # local minima, gradient descent, pointing from white to black (1->0)
#         ly[ii] = -1.*scalespace.get('ly')
#
#     feature_ori = []
#     for ii, by in enumerate(argmaxgrad[1]):
#         scale_index = np.where(scale_range == tnew[ii])
#         bx = argmaxgrad[0][ii]
#         sigma_ori = sigma_ori_times*np.sqrt(tnew[ii])  # 1.5 is gaussian sigma for orientation assignment (lambda_ori
#         #  in Re13)
#         radius = int(np.floor(window_ori_radtimes*sigma_ori))  # radius of the (squared) region used in orientation
#         # assignment
#         ori = orientation_hist(l[scale_index[0]], bx, by, lx[scale_index[0]], ly[scale_index[0]],
#                                radius, sigma_ori, n_bins_ori, peak_ratio, smooth_cycles, image_original=image)
#         feature_ori.append(ori)
#
#         plt.figure()
#         plt.imshow(image.T, interpolation='none', cmap='gray', origin='low')
#         plt.hold(True)
#         for i in ori:
#             plt.arrow(bx, by, 10*np.cos(i), 10*np.sin(i), head_width=1, head_length=1, fc='r', ec='r')
#             plt.show()
#
#     return feature_ori


def orientation_hist(image, bx, by, lx, ly, radius, sigma_ori, n_bins_ori, peak_ratio=0.8, smooth_cycles=2,
                     image_original=None, plot_graphics=True):
    """
    Computes a gradient orientation histogram at a specified pixel
    compute gradient values, orientations and the weights over the pixel neighborhood
    """

    from matplotlib.colors import LogNorm

    hist = np.zeros(shape=n_bins_ori)
    if bx-np.sqrt(2)*radius < 0 or bx+np.sqrt(2)*radius > image.shape[1] or by-np.sqrt(2)*radius < 0 or \
                            by+np.sqrt(2)*radius > image.shape[2]:
        return []
    grad_x = lx[bx-radius:bx+radius+1, by-radius:by+radius+1]
    grad_y = ly[bx-radius:bx+radius+1, by-radius:by+radius+1]
    locgrid_x, locgrid_y = np.meshgrid(np.linspace(bx - radius, bx + radius, 2*radius+1),
                                       np.linspace(by - radius, by + radius, 2*radius + 1))
    weights = np.exp(-((locgrid_x-bx)**2 + (locgrid_y-by)**2)/(2*sigma_ori**2))

    ori = np.arctan2(grad_x, grad_y)
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
        plt.imshow(image_original.T, interpolation = 'none', cmap = 'gray', origin = 'lower', norm = LogNorm())
        plt.hold(True)
        for i in range(0, 2*radius+1,4):
            y = by + i - radius
            for j in range(0, 2*radius+1,4):
                x = bx + j - radius  # col
                ax.arrow(x, y, weights[j, i]*grad_y[j, i], weights[j, i]*grad_x[j, i],
                         head_width=0.2, head_length=0.2, fc='r', ec='r')
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
        return []

    # image = image_test
    # plt.close('all')
    # n_hist = 16
    # n_bins_descr = 8
    # window_descr_radtimes = 4
    # argmaxgrad = feature.get('argmaxgrad')
    # tnew = feature.get('tnew')  # 1d-array
    # scale_range = feature.get('scale_range')
    # orientation = feature.get('orientation')
    # sigma_descr_times = 1.5
    # sigma_descr = sigma_descr_times * np.sqrt(tnew[0])
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
    grid_x, grid_y = np.meshgrid(np.arange(0, image.shape[0]), np.arange(0, image.shape[1]))
    if bx - np.sqrt(2)*radius < 0 or bx + np.sqrt(2)*radius > image.shape[0] or by - np.sqrt(2)*radius < 0 or by + \
            np.sqrt(2)*radius > image.shape[1]:
        return []

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
        #print(pos_ref[0][63])
        grad_x = lx[pos[1], pos[0]]; grad_y = ly[pos[1], pos[0]]
        ori_all = (np.arctan2(grad_x, grad_y) - ori)*180.//np.pi % 360

        weights = np.exp(-((pos[0]-by)**2 + (pos[1]-bx)**2)/(2*radius**2))
        contr = np.multiply(weights, np.sqrt(grad_x**2 + grad_y**2))  # contribution to the histogram

        # build the histogram
        hist = np.zeros(shape=(n_hist*n_bins_descr))
        ori_all_bins = np.floor(ori_all / (360 / n_bins_descr))
        # print '\n******************ori=', ori
        # print np.min(hist_ind), np.min(np.arctan2(grad_x, grad_y) - ori), np.min(ori_all), np.min(ori_all_bins)
        # print np.max(hist_ind), np.max(np.arctan2(grad_x, grad_y) - ori), np.max(ori_all), np.max(ori_all_bins)
        # print 'hist_ind[0]=', hist_ind[0]
        # print ori_all[0], ori_all_bins[0]
        bins_temp = (n_bins_descr * hist_ind + ori_all_bins).astype(int)
        for i, bins in enumerate(bins_temp):
            #print i, bins
            hist[bins] += contr[i]

        # normalize descriptor vectors
        hist /= np.linalg.norm(hist)
        hist[np.where(hist > threshold_sat)] = threshold_sat  # we reduce the influence of large gradient
                # magnitudes by thresholding the values in the unit feature vector to each be no larger than
                # 0.2, and then renormalizing to unit length.
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

            plt.figure(); plt.title('center = (%d,%d), ori = %f' % (bx,by,ori))
            ax = plt.subplot(1, 2, 1)
            plt.imshow(image.T, interpolation='none', cmap='gray', origin='lower', norm=LogNorm())
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
                plt.bar(np.arange(n_bins_descr), hist[n_bins_descr*ind:n_bins_descr*ind+n_bins_descr], color='0.5',
                        align='center')
                plt.xticks(np.linspace(0, n_bins_descr, 5), ['0', '\pi/2', '\pi', '3\pi/2', '2\pi'])
                plt.ylim([0, max_ori])
                plt.hold(True)


    return histogram
