import os, matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams["text.usetex"] = True; matplotlib.rcParams['font.family'] = 'serif'  # configure latex plots


def _adjust_bounds(ax, points):

    ptp_bound = points.ptp(axis=0)
    ax.set_xlim(points[:, 0].min() - 0.1 * ptp_bound[0],
                points[:, 0].max() + 0.1 * ptp_bound[0])
    ax.set_ylim(points[:, 1].min() - 0.1 * ptp_bound[1],
                points[:, 1].max() + 0.1 * ptp_bound[1])


def voronoi_plot_2d(vor, **kw):
    """
    Plot the given Voronoi diagram in 2-D, based on scipy.spatial.voronoi_plot_2d

    Input:
    ----------
    vor : scipy.spatial.Voronoi instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on
    show_points: bool, optional
        Add the Voronoi points to the plot.
    show_vertices : bool, optional
        Add the Voronoi vertices to the plot.
    line_colors : string, optional
        Specifies the line color for polygon boundaries
    line_width : float, optional
        Specifies the line width for polygon boundaries
    line_alpha: float, optional
        Specifies the line alpha for polygon boundaries

    Output:
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot

    """
    from matplotlib.collections import LineCollection
    from distutils.util import strtobool

    fig, ax = plt.subplots()

    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    if strtobool(kw.get('show_vertices', True)):
        ax.plot(vor.vertices[:, 0], vor.vertices[:, 1], 'o')

    line_colors = kw.get('line_colors', 'k')
    line_width = kw.get('line_width', 1.0)
    line_alpha = kw.get('line_alpha', 1.0)

    line_segments = []
    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):  # exclude lines to infinity
            line_segments.append([(x, y) for x, y in vor.vertices[simplex]])

    lc = LineCollection(line_segments,
                        colors=line_colors,
                        lw=line_width,
                        linestyle='solid')

    lc.set_alpha(line_alpha)
    ax.add_collection(lc)
    ptp_bound = vor.points.ptp(axis=0)

    if strtobool(kw.get('show_unbounded_cells', True)):
        line_segments = []
        center = vor.points.mean(axis=0)
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            if np.any(simplex < 0):
                i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

                t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[pointidx].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[i] + direction * ptp_bound.max()

                line_segments.append([(vor.vertices[i, 0], vor.vertices[i, 1]),
                                      (far_point[0], far_point[1])])

        lc = LineCollection(line_segments,
                            colors=line_colors,
                            lw=line_width,
                            linestyle='dashed')
        lc.set_alpha(line_alpha)
        ax.add_collection(lc)
    _adjust_bounds(ax, vor.points)

    return ax.figure


def plot_points(points, title='', plot_axis='on'):

    fig, ax = plt.subplots()
    ax.plot(points[:, 0], points[:, 1], 'k.', markersize=0.5, alpha=0.5)
    fig.show(); ax.hold(True)

    plt.title(title)
    ax.set_xlim(points[:, 0].min(), points[:, 0].max()); ax.set_ylim(points[:, 1].min(), points[:, 1].max())
    ax.set_aspect('equal', adjustable='box')

    if plot_axis is 'off':
        fig.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    ax.hold(False)


def compute_parameters(vor, dict_inputfile):
    """
    Compute quatitative paramters form the voronoi diagram. Modifies the input object class vor with new attributes

    Input:
    ----------
    vor : Voronoi object class
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    original_pixel_size = dict_inputfile.get('original_pixel_size')
    compute_areas(vor, original_pixel_size)
    compute_densities(vor)


def compute_densities(vor):
    """
    Compute densities at the 0-th rank (see Levet et al., 2015). First run vprocessing.compute_areas

    Input:
    ----------
    vor: Voronoi object class

    Output:
    ----------
    new attribut vor.densities_zero and vor.densities_average
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")
    if not hasattr(vor, 'areas'):
        raise ValueError("Requires vor.areas attribute! Please, run vprocessing.compute_areas(vor).")

    vor.densities_zero = 1./vor.areas
    vor.densities_average = vor.points.shape[0]/vor.areas_total


def compute_areas(vor, original_pixel_size):
    """
    Compute voronoi areas (vareas) given  infinite voronoi regions in a 2D diagram to finite
    regions. The area is at the 0-th rank (see Levet et al., 2015). Modifies the input object class vor with new
    attribute areas (vor.areas = 1d-array corresponding to each  vor.point_region. -1 indicates unbounded region or
    finite region outside ROI)

    Input
    ----------
    vor : Voronoi object class

    Output:
    ----------
    new attribut vor.areas and vor.areas_total
    """
    from shapely.geometry import MultiPoint, Point, Polygon

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    areas = np.empty(shape=vor.point_region.shape)
    mask = MultiPoint([Point(i) for i in vor.points]).convex_hull.buffer(0)
    for p1, region in enumerate(vor.point_region):
        vertices_label = vor.regions[region]  # labels of the voronoi vertex forming the boundary of that region,
        # containing point p1
        if -1 in vertices_label:  # vertex label -1 means infinite region
            areas[p1] = float('inf')  # None #-1  # None area -1 means unbounded region
        else:
            polygon_vertices = vor.vertices[vertices_label]
            poly = Polygon(polygon_vertices)
            intersec = poly.intersection(mask)
            poly_x, poly_y = poly.exterior.xy
            intersec_x, intersec_y = intersec.exterior.xy
            areas[p1] = 0.5 * np.abs(np.dot(polygon_vertices[:, 0], np.roll(polygon_vertices[:, 1], -1))
                                     - np.dot(polygon_vertices[:, 1], np.roll(polygon_vertices[:, 0], -1)))
            if set(list(poly_x)) == set(list(intersec_x)) and set(list(poly_y)) == set(list(intersec_y)):
                # Shoelace formula
                areas[p1] = 0.5*np.abs(np.dot(polygon_vertices[:, 0], np.roll(polygon_vertices[:, 1], -1))
                                       - np.dot(polygon_vertices[:, 1], np.roll(polygon_vertices[:, 0], -1)))
            else:  # polygon intersection convex hull is not polygon
                areas[p1] = float('inf')  # -1 area -1 means voronoi region outside the convex hull (~ROI)

    vor.areas = original_pixel_size**2*areas   # in nm2
    vor.areas_total = float(np.sum(vor.areas[(vor.areas > 0) & (vor.areas < float('inf'))]))  # in nm2


def plot_areas(vor, thr=None, plot_axis='on', show_points=True, hold=False):
    """
    This function plots:
    1) the areas of the bounded voronoi regions within ROI, contained in vor.areas.
    2) if threshold (thr) is given, it also plots voronoi pologons with area smaller than threshold.

    Input:
    -----------------------------
    thr (float): optional maximum area (in nm2) to generate a thresholded voronoi

    """
    from matplotlib.collections import PolyCollection

    if not hasattr(vor, 'areas'):
        raise ValueError("Requires vor.areas attribute! Please, run vprocessing.compute_areas(vor).")

    max_area = np.max(vor.areas[vor.areas < float('inf')])
    polygons, polygons_thresholded = [], []
    p = 0
    for p1, area in enumerate(vor.areas):
        if (area >= 0) and (area < float('inf')):
            vertices_label = vor.regions[vor.point_region[p1]]
            polygons.append([(x, y) for x, y in vor.vertices[vertices_label]])
            p += 1
            if area < thr:
                polygons_thresholded.append(p-1)
    fig, ax = plt.subplots()

    lsc = 1 / 1e-04    # Jerome's colour:
    sc = 999 / np.log10(max_area * lsc)
    colour = (1000-sc*np.log10(lsc*vor.areas[(vor.areas >= 0) & (vor.areas < float('inf'))]))
    # colour = 1-np.log10(vor.areas[vor.areas >= 0])   # my colour:
    coll = PolyCollection(polygons, array=colour, edgecolors='none'); ax.add_collection(coll), ax.autoscale_view()
    fig.colorbar(coll, ax=ax)   # Add a colorbar for the PolyCollection
    plt.hold(True)
    if show_points is True:
        ax.plot(vor.points[:, 0], vor.points[:, 1], 'k.', markersize=2)
    if plot_axis is 'off':
        plt.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    ax.set_xlim(vor.points[:, 0].min(), vor.points[:, 0].max())
    ax.set_ylim(vor.points[:, 1].min(), vor.points[:, 1].max())
    ax.set_aspect('equal', adjustable='box')
    fig.hold(hold)

    if thr is not None:
        fig_th, ax_th = plt.subplots()
        polygons_th = [polygons[i] for i in polygons_thresholded]
        import matplotlib.colors as colors
        coll_th = PolyCollection(polygons_th, color='grey', edgecolors='none')
        ax_th.add_collection(coll_th)
        ax_th.set_xlim(vor.points[:, 0].min(), vor.points[:, 0].max())
        ax_th.set_ylim(vor.points[:, 1].min(), vor.points[:, 1].max())
        ax_th.set_aspect('equal', adjustable='box')

    return fig, ax  # .figure


def plot_densities(vor, thr=None, plot_axis='on', show_points=True, cmap='jet', norm='linear', hold=False):
    """
    This function plots:
    1) the zero-rank densities of the bounded voronoi regions within ROI, contained in
    vor.densities_zero
    2) if threshold is given, it also plots voronoi pologons with densities larger than threshold.

    Input:
    -----------------------------
    threshold (float): optional minimum density (in nm-2) to generate a thresholded voronoi

    """

    from matplotlib.collections import PolyCollection
    from matplotlib import cm, colors

    if not hasattr(vor, 'densities_zero'):
        raise ValueError("Requires vor.densities_zero attribute! Please, run vprocessing.compute_densities(vor).")

    polygons, polygons_thresholded = [], []
    p = 0
    for p1, density in enumerate(vor.densities_zero):
        if density > 0:
            vertices_label = vor.regions[vor.point_region[p1]]
            polygons.append([(x, y) for x, y in vor.vertices[vertices_label]])
            p += 1
            if density >= thr:
                polygons_thresholded.append(p-1)
    fig, ax = plt.subplots()

    if cmap is 'jet': cmap = plt.cm.jet
    if cmap is 'gray': cmap = plt.cm.gray

    if norm is 'linear':
        values = vor.densities_zero  # vor.densities_zero[vor.densities_zero > 0]
    if norm is 'log':
        values = np.log10(vor.densities_zero[vor.densities_zero > 0])
    cnorm = colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    scalarmap = cm.ScalarMappable(norm=cnorm, cmap=cmap); scalarmap.set_array(values)
    colour = [scalarmap.to_rgba(ii) for ii in values]
    coll = PolyCollection(polygons, edgecolors='none'); ax.add_collection(coll), ax.autoscale_view()
    coll.set_facecolors(colour)
    # cbar = fig.colorbar(coll, ax=ax)   # Add a colorbar for the PolyCollection
    # cbar.ax.set_ylabel('zero-rank density [nm$^{-2}$]', rotation=270); cbar.ax.set_xlabel('$log_{10}$')
    plt.hold(True)
    if show_points is True:
        ax.plot(vor.points[:, 0], vor.points[:, 1], 'k.', markersize=0.5)

    if plot_axis is 'off':
        plt.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    ax.set_xlim(vor.points[:, 0].min(), vor.points[:, 0].max())
    ax.set_ylim(vor.points[:, 1].min(), vor.points[:, 1].max())
    ax.set_aspect('equal', adjustable='box')
    fig.hold(hold)

    if thr is not None:
        fig_th, ax_th = plt.subplots()
        polygons_th = [polygons[i] for i in polygons_thresholded]
        coll_th = PolyCollection(polygons_th, color='grey', edgecolors='none')
        ax_th.add_collection(coll_th)  # , ax_th.autoscale_view()
        ax_th.set_xlim(vor.points[:, 0].min(), vor.points[:, 0].max())
        ax_th.set_ylim(vor.points[:, 1].min(), vor.points[:, 1].max())
        ax_th.set_aspect('equal', adjustable='box')

    return fig, ax  # .figure


def threshold(vor, thr=None):
    """
    This function thresholds the point pattern based on the voronoi densities:

    Input:
    -----------------------------
    threshold (float): optional minimum density (in nm-2) to generate a thresholded voronoi

    Outpu:
    ----------------------------
    new attribute vor.points_thresholded = vor.points with densities larger than threshold

    """
    if not hasattr(vor, 'densities_zero'):
        raise ValueError("Requires vor.densities_zero attribute! Please, run vprocessing.compute_densities(vor).")

    points_thresholded = []
    for p1, density in enumerate(vor.densities_zero):  # p1 is region index, no point
        if density > 0 and density >= thr:
            points_thresholded.append(vor.points[p1])

    vor.points_thresholded = np.asarray(points_thresholded)

    return vor


def densities_interpolate(vor, scale_pixel_size, interpolate_method='nearest', fill_value=0.0,
                          density_transform='linear'):
    """
    Given a vornoi-based zero-rank density map, the natural interpolation (Sibson, 1980) on a regular grid with
    analysis_pixel_size grid size is performed. See Andronov et al., 2016 and Ledoux and Gold, 2005.

    Output:
    -----------------------
    densities_image (2d array)
    """
    from scipy.interpolate import griddata
    from matplotlib.colors import LogNorm

    grid_x, grid_y = np.mgrid[np.min(vor.points[:, 0]):np.max(vor.points[:, 0]):scale_pixel_size,
                              np.min(vor.points[:, 1]):np.max(vor.points[:, 1]):scale_pixel_size]
    densities_image = griddata(vor.points, vor.densities_zero, (grid_x, grid_y),
                               method=interpolate_method, fill_value=fill_value)
    # fig, ax = plt.subplots()
    # plt.imshow(densities_image.T, extent=(np.min(grid_x[:, 0]), np.max(grid_x[:, 0]), np.min(grid_y[0, :]),
    #            np.max(grid_y[0, :])), interpolation='none', cmap='jet', origin='lower', norm=LogNorm())
    # cbar = plt.colorbar()
    # cbar.ax.set_ylabel('zero-rank density [nm$^{-2}$]', rotation=270); cbar.ax.set_xlabel('$log_{10}$')
    # plt.show()

    if density_transform is not 'linear':
        if np.any(densities_image <= 0):
            densities_image[np.where(densities_image <= 0)] = None
        densities_image = np.log10(densities_image)

    return densities_image


def plot_feature(vor, feature, dict_sift, plot_axis='on', show_points=False, cmap='jet', norm=None,
                 blob_color='strength', ori_color=None, ori_cmap=None):
    """
    Function for plotting scale-space feature on Vornoi tessellation (STORM pixel size)

    Input:
    --------
    feature (dictionary)
    cmap: map of the image density map (if not just points)
    norm (string): 'log', 'lin'. Scale of imshow. E.g.,
        if norm is 'log':  plt.imshow(image.T, interpolation='none', cmap='jet', origin='lower', norm=LogNorm())
    blob_color = 'strength', 'scale', 'class'
    ori_color (list) = if not None, for all orientations, assigned classes
    ori_cmap = if not None, discrete cmap for vocabulary - output clustering/unsupervised learing.

    """
    import matplotlib.colors as colors
    from collections import Counter

    feature_name = dict_sift.get('feature_name')
    scale_pixel_size = dict_sift.get('scale_pixel_size')
    argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    scale = feature.get('scale')  # 1d-array
    strength = feature.get('strength')
    orientation = feature.get('orientation', [])

    # fig, ax = plot_densities(vor, plot_axis=plot_axis, show_points=show_points, cmap=cmap, norm=norm, hold=True)
    fig, ax = plt.subplots()
    ax.plot(vor.points[:, 0], vor.points[:, 1], 'k.', markersize=0.5, alpha=0.5); ax.hold(True)
    ax.set_xlim(vor.points[:, 0].min(), vor.points[:, 0].max()); ax.set_ylim(vor.points[:, 1].min(), vor.points[:, 1].max())
    ax.set_aspect('equal', adjustable='box'); ax.hold('True')

    ox = np.min(vor.points[:, 0]); oy = np.min(vor.points[:, 1])

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
                    plt.plot([ox + bx * scale_pixel_size, ox + bx * scale_pixel_size +
                              np.sqrt(scale[ii]) * 1 * 1.5 * np.cos(ori) * scale_pixel_size],
                             [oy + by * scale_pixel_size, oy + by * scale_pixel_size +
                              np.sqrt(scale[ii]) * 1 * 1.5 * np.sin(ori) * scale_pixel_size],
                              color=o_color, linewidth=0.5)
                if len(orientation[ii]) > 0 and ori_color is not None: mean = Counter(mean).most_common(1)[0][0]
            # # plot blobs - detected features
            if blob_color == 'strength':
                b_color = blob_cmap.to_rgba(strength[ii])
            elif blob_color == 'scale':
                b_color = blob_cmap.to_rgba(np.sqrt(scale[ii]) * 2 * 1.5)
            elif blob_color == 'class':
                if len(orientation[ii]) == 0:
                    b_color = 'None'
                else:
                    b_color = blob_cmap(mean)
            ax = plt.plot(ox + (ucirc[0, :] * np.sqrt(scale[ii]) * 1*1.5 + bx)*scale_pixel_size,
                          oy + (ucirc[1, :] * np.sqrt(scale[ii]) * 1*1.5 + by)*scale_pixel_size,
                          color=b_color, linewidth=0.5)
        # fig.colorbar(blob_cmap, label='max$_t\{\Delta_{\gamma-norm}\}$')
        # fig.colorbar(scalarmap, label='3$\sqrt{t}$ [pixel - original]')

    if plot_axis is 'off':
        plt.axis(plot_axis)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
    plt.hold(False)

    return ax


def localizations_feature(vor, feature, dict_sift):
    """
    Count number of localizations within intensity-dependent and feature-based clusters (default: blobs)

    Input:
    -------------

    Output:
    -------------
    number_localizations (1d-array): following the order of argmaxgrad, each position contains the number of
    localizations within that feature (e.g., blob with radius 2*\sqrt{t})
    """

    scale_pixel_size = dict_sift.get('scale_pixel_size')
    argmaxgrad = feature.get('argmaxgrad')  # tuple of (argmaxgrad[0], argmaxgrad[1]) = (ndarray, ndarray) = (col, row)
    scale = feature.get('scale')  # 1d-array

    ox = np.min(vor.points[:, 0]); oy = np.min(vor.points[:, 1])
    number_localizations = np.zeros(shape=(len(argmaxgrad[0]),))
    for ii, blob_y in enumerate(argmaxgrad[1]):
        blob_x = argmaxgrad[0][ii]
        radius = (np.sqrt(scale[ii]) * 1.5)*scale_pixel_size
        count = np.where(np.sqrt((vor.points[:, 0] - (ox + blob_x*scale_pixel_size))**2 +
                         (vor.points[:, 1] - (oy + blob_y*scale_pixel_size))**2) <= radius)
        number_localizations[ii] = len(count[0])

    return number_localizations
