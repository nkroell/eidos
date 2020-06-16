# -*- coding: utf-8 -*-
"""

@author: Nils Kroell
"""

import numpy as np
import cv2
from skimage import measure, morphology, transform
from scipy import spatial, stats
import matplotlib.pyplot as plt
import pandas as pd


def extract_all_shape_measurements_2d(bw, dalpha=9, return_statistical_lengths_distributions=False, return_all_chords=False):
    """ Calculates 2D shape measurements of particles in a binary image.
    
        Args:
            bw (ndarry, bool): Binary image.
            dalpha (int/float): Angle in degree [°] to turn bw each iteration. Should be a fraction of 180° (18ß°/n), where n
                is a natural number.
            return_statistical_lengths_distributions (boolean): If True then the statistical length are returned.
                Default is False.
            return_all_chords (boolean): If True then all found chords are returned.
                Default is False.
        
        Returns:
            df_shape_measurements (Pandas Dataframe): A pandas dataframe with all extracted shape measurements (columns) 
                for all found particles (rows).
            dfs_statistical_lengths_distributions (list of Pandas Dataframes): A list of n pandas dataframes, where each
                panda dataframe includes the statistical lengths (columns) for all iterated angles (0°:dalpha:180°).
                Only returned if return_statistical_lengths_distributions is True.
            return_all_chords (ndarray): A numpy array of with all found chords for all angles.
                Only returned if return_all_chords is True.
    """
    
    # Initialize results
    all_chords = np.array([])
    dfs_shape_measurements = []
    dfs_statistical_lengths_distributions = []
    
    # Label binary image to extract single objects
    labels, n_objects = measure.label(bw, return_num=True)    
    # Iterate over all objects
    for i in range(1,n_objects+1):
        bw_single = np.zeros(bw.shape)
        bw_single[labels == i] = True
        df_shape_measurements_single_particle, df_statistical_lengths_distributions_single_particle, all_chords_single_particle = _shape_measurements_2d_single_object(bw_single, dalpha)
        dfs_shape_measurements.append(df_shape_measurements_single_particle)
        dfs_statistical_lengths_distributions.append(df_statistical_lengths_distributions_single_particle)
        all_chords = np.append(all_chords, all_chords_single_particle)

    # Merge dataframes
    if len(dfs_shape_measurements) > 1:
        df_shape_measurements = pd.concat(dfs_shape_measurements, ignore_index=True)
    else:
        df_shape_measurements = dfs_shape_measurements[0]

        
    
    
    # Return depending on settings
    if return_statistical_lengths_distributions == False and return_all_chords == False:
        return df_shape_measurements
    elif return_statistical_lengths_distributions == True and return_all_chords == False:
        return df_shape_measurements, dfs_statistical_lengths_distributions
    elif return_statistical_lengths_distributions == False and return_all_chords == True:
        return df_shape_measurements, all_chords
    else:
        return df_shape_measurements, dfs_statistical_lengths_distributions, all_chords


def _shape_measurements_2d_single_object(bw, dalpha):
    """ Calculates 2D rotation and translation invariant shape measurements of a binary image, containing one object.
    
    Args:
        bw (ndarray, bool): Binary image.
        
    Returns:
    
        df_shape_measurements_single_particle (Pandas Dataframe): Pandas dataframe with one row with shape measurements.
        df_statistical_lengths_distributions_single_particle (Pandas Dataframe): Pandas dataframe with statistical length at their
            measured angle.
        all_chords_single_particle (ndarray): Array of all found choords (for all rotations).
    """
   
    # assumes: only one particle per binary image

    ##### MARCRO DESCRIPTORS ##### 
    perimeter, area, filled_area, convex_area, major_axis_length, minor_axis_length, centroid, coords, bw_cropped , bw_convex_cropped = calc_skimage_measurements(bw)

    # Contour (for further processing)
    contour = calc_contour_list(bw)

    # Circles
    max_inclosing_circle_center, max_inclosing_circle_radius = calc_max_inclosing_circle(bw)
    min_enclosing_circle_center, min_enclosing_circle_radius = calc_min_enclosing_circle(contour)
    circumscribing_circle_radius, inscribing_circle_radius = calc_circumscribing_and_inscribing_circle(centroid, contour)
    
    # use diameters instead of radii:
    max_inclosing_circle_diameter = 2 * max_inclosing_circle_radius
    min_enclosing_circle_diameter = 2 * min_enclosing_circle_radius
    circumscribing_circle_diameter = 2 * circumscribing_circle_radius
    inscribing_circle_diameter = 2 * inscribing_circle_radius

    # Statistical length
    # distributions
    feret_diameters, martin_diameters, nassenstein_diameters, max_chords, all_chords, measured_angles  = calc_statistical_length_distribution(bw_cropped, daplha=dalpha)
    # distribution parameters
    max_feret, min_feret, median_feret, mean_feret, mode_feret, std_feret = calc_distribution_parameters(feret_diameters)
    max_martin, min_martin, median_martin, mean_martin, mode_martin, std_martin = calc_distribution_parameters(martin_diameters)
    max_nassenstein, min_nassenstein, median_nassenstein, mean_nassenstein, mode_nassenstein, std_nassenstein = calc_distribution_parameters(nassenstein_diameters)
    max_max_chords, min_max_chords, median_max_chords, mean_max_chords, mode_max_chords, std_max_chords = calc_distribution_parameters(max_chords)
    max_all_chords, min_all_chords, median_all_chords, mean_all_chords, mode_all_chords, std_all_chords = calc_distribution_parameters(all_chords)

    # Main and maximum dimensions
    x_max, y_max = max_dimensions(max_chords, measured_angles)
    width_min_bb, height_min_bb, center_bb, cornerpoints_min_bb = calc_min_bounding_box(coords)

    # Equal diameters
    area_equal_diameter = calc_area_equal_diameter(area)
    perimeter_equal_diameter = calc_perimeter_equal_diameter(perimeter)

    # Geodatic length and thickness
    geodeticlength, thickness = calc_geodeticlength_and_thickness(area, perimeter)


    ##### MESODESCRIPTORS ##### 
    convex_perimeter = calc_convex_perimeter(bw_convex_cropped)

    # Measurements based on erosion
    n_erosions_binary_image = calc_n_erosions_to_erase_binary_img(bw_cropped)
    n_erosions_complement = calc_n_erosions_to_erase_binary_complement(bw_cropped, bw_convex_cropped)


    ##### MICRODESCRIPTORS ##### 
    fractal_dimension_boxcounting_method = calc_fractal_dimension_boxcounting_method(bw_cropped)
    fractal_dimension_perimeter_method = calc_fractal_dimension_perimeter_method(contour, max_feret)

    shape_measurements = {
        "perimeter": perimeter,
        "convex_perimeter": convex_perimeter,
        "area": area,
        "filled_area": filled_area,
        "convex_area": convex_area,
        "major_axis_length": major_axis_length,
        "minor_axis_length": minor_axis_length,
        "max_inclosing_circle_diameter": max_inclosing_circle_diameter,
        "min_enclosing_circle_diameter": min_enclosing_circle_diameter,
        "circumscribing_circle_diameter": circumscribing_circle_diameter,
        "inscribing_circle_diameter": inscribing_circle_diameter,
        "x_max": x_max,
        "y_max": y_max,
        "width_min_bb": width_min_bb,
        "height_min_bb": height_min_bb,
        "area_equal_diameter": area_equal_diameter, 
        "perimeter_equal_diameter": perimeter_equal_diameter,
        "geodeticlength": geodeticlength,
        "thickness": thickness,
        "n_erosions_binary_image": n_erosions_binary_image,
        "n_erosions_complement": n_erosions_complement,
        "fractal_dimension_boxcounting_method": fractal_dimension_boxcounting_method,
        "fractal_dimension_perimeter_method": fractal_dimension_perimeter_method,    
        "max_feret": max_feret,
        "min_feret": min_feret,
        "median_feret": median_feret,
        "mean_feret": mean_feret,
        "mode_feret": mode_feret,
        "std_feret": std_feret,
        "max_martin": max_martin,
        "min_martin": min_martin,
        "median_martin": median_martin,
        "mean_martin": mean_martin,
        "mode_martin": mode_martin,
        "std_martin": std_martin,
        "max_nassenstein": max_nassenstein,
        "min_nassenstein": min_nassenstein,
        "median_nassenstein": median_nassenstein,
        "mean_nassenstein": mean_nassenstein,
        "mode_nassenstein": mode_nassenstein,
        "std_nassenstein": std_nassenstein,
        "max_max_chords": max_max_chords,
        "min_max_chords": min_max_chords,
        "median_max_chords": median_max_chords,
        "mean_max_chords": mean_max_chords,
        "mode_max_chords": mode_max_chords,
        "std_max_chords": std_max_chords,
        "max_all_chords": max_all_chords,
        "min_all_chords": min_all_chords,
        "median_all_chords": median_all_chords,
        "mean_all_chords": mean_all_chords,
        "mode_all_chords": mode_all_chords,
        "std_all_chords": std_all_chords 
    }

    statistical_lengths_distributions = {
        "measured_angle": measured_angles,
        "feret_diameter": feret_diameters,
        "martin_diameter": martin_diameters,
        "nassenstein_diameter": nassenstein_diameters,
        "max_chord": max_chords
    }

    df_shape_measurements_single_particle = pd.DataFrame.from_dict(shape_measurements)
    df_statistical_lengths_distributions_single_particle = pd.DataFrame.from_dict(statistical_lengths_distributions)
    
    all_chords_single_particle = all_chords

    return df_shape_measurements_single_particle, df_statistical_lengths_distributions_single_particle, all_chords_single_particle


def calc_skimage_measurements(bw):
    """ Calculates particle shape measurements based on skimage.measure.regionprops.
    
    Args:
        bw (ndarray, boolean): Binary image describing the particle shape. (True = particle, False = background)
        
    Returns:
        area (float): The area of the particle
        perimeter (float): The perimeter of the particle.
        convex_area (float): The convex area of the particle.
        major_axis_length (float): The major axis length of the particle.
        minor_axis_length (float): The minor axis length of the particle.
        centroid (ndarray): The centroid of the particle in the form [x, y].
        coords (ndarray): A list of coordinates of the particles (indexes of binary images)
            in the form [[x0, y0], [x1, y1], ...].
    """
    
    labels = measure.label(bw)
    props = measure.regionprops_table(labels, properties=('area', 'filled_area', 'convex_area', 'major_axis_length', 'minor_axis_length', 'perimeter',
                                                          'coords', 'centroid', 'convex_image', 'image'))
    
    perimeter = props['perimeter'][0]
    
    area = props['area'][0]
    filled_area = props['filled_area'][0]
    convex_area = props['convex_area'][0]
    
    major_axis_length = props['major_axis_length'][0]
    minor_axis_length = props['minor_axis_length'][0] 
    
    centroid = np.zeros((2,))
    centroid[0], centroid[1] = props['centroid-0'][0], props['centroid-1'][0]
    
    coords = props['coords'][0]
    
    bw_cropped = props['image'][0]
    bw_convex_cropped = props['convex_image'][0]
    
    return perimeter, area, filled_area, convex_area, major_axis_length, minor_axis_length, centroid, coords, bw_cropped , bw_convex_cropped 


def calc_contour_list(bw):
    """ Returns a list of **external** contour points of a particle shape.
    
    Args:
        bw (ndarray, boolean): A binary image of the particle shape.
        
    Returns:
        contour (ndarray, int): A list of indexes of the contourpoints of the particle
            in the form [[x0, y0], [x1, y1], [x2, y2]]    
    """
    
    # find contour
    contour_cv2, hierarchy = cv2.findContours(bw.astype('uint8'), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    
    # transform list of contour points and change coordinate system to (x,y)
    contour_cv2 = contour_cv2[0].reshape((-1,2))
    contour = np.zeros(contour_cv2.shape)
    contour[:,0], contour[:,1] = contour_cv2[:,1], contour_cv2[:,0]
    
    return contour


def calc_max_inclosing_circle(bw):
    """ Calculates the maximum inclosing circle of a particle shape.
    
    Args:
        bw (ndarray, boolean): binary image of the particle (true = particle, false = background).
        
    Returns:
        max_inclosing_circle_center (ndarray, float): center position of the maximum inclosing circle in form [x, y].
        max_inclosing_circle_radius (float): Radius of the maximum inclosing circle . 
    """
    distance_map = cv2.distanceTransform(bw.astype('uint8'), cv2.DIST_L2, cv2.DIST_MASK_3)
    radius= np.amax(distance_map)
    idx_max_distance = np.where(distance_map==radius)
    center = np.zeros((2,))
    center[0], center[1] = idx_max_distance[0][0], idx_max_distance[1][0]
    
    return center, radius


def calc_min_enclosing_circle(contour):
    """ Calculates the minimum enclosing circle of a particle based on its list of coordinates
    
    Args:
        contour (ndarray, int/float): List of contour points describing the particle in the form [[x0, y0], [x1, y1], ...].
              Note: If float => converted to int.
    Returns:
        center (ndarray, float): The center of the minium enclosing circle.
        radius (float): The radius of the minium enclosing circle.    
    """
    
    # transform coordinates to cv2 representation:
    coords_cv2 = np.zeros(contour.shape)
    coords_cv2[:,0] = contour[:,1]
    coords_cv2[:,1] = contour[:,0]
    
    coords_cv2 = (coords_cv2.reshape((-1,2,1))).astype('int32')
    
    # find min enclosing circle
    center_cv2, radius = cv2.minEnclosingCircle(coords_cv2)
    
    # convert to numpy array and transform to (x,y) coordinate system
    center = np.zeros((2,))
    center[0], center[1] = center_cv2[1], center_cv2[0]
    
    return center, radius


def calc_circumscribing_and_inscribing_circle(centroid, contour):
    """ Calculates the circumscribing and inscribing circle of a particle shape.
    
    Note: Uses scipy.spatial.distance.cdist to calculate distances.
    
    Args:
        centroid (ndarray): Centroid of the particle in the form [x,y].
        contour (ndarray): Array describing the contour of the particle. Contour must have shape (n_contourpoints,2) [[x0, y0], [x1, y1], ...]
    
    Returns:
        circumscribing_circle_radius (float): The radius of the circumscribing circle.
            The circumscribing circle touches the contour from the outside, has the same centroid as the contour and minimum area.
        inscribing_circle_radius (float): The radius of the inscribing circle.
            The circumscribing circle touches the contour from the inside, has the same centroid as the contour and maximum area.
    """
    
    distances = spatial.distance.cdist(centroid.reshape((-1,2)), contour, 'euclidean')
    distances = distances.reshape((-1))
    circumscribing_circle_radius = max(distances)
    inscribing_circle_radius = min(distances)
    
    return circumscribing_circle_radius, inscribing_circle_radius


def calc_area_equal_diameter(area):
    """ Calculates the areaequal diameter of a particle
    Args:
        area (float): Area of the particle.
        
    Returns:
        area_equal_diameter (float): Areaequal diameter of the particle.
    """
    
    return np.sqrt(4*area/np.pi)


def calc_perimeter_equal_diameter(perimeter):
    """ Calculates the perimeter equal diameter of a particle
    Args:
        perimeter (float): Perimeter of the particle.
        
    Returns:
        perimeter_equal_diameter (float): Perimeterequal diameter of the particle.
    """
    return perimeter/np.pi


def calc_geodeticlength_and_thickness(area, perimeter):
    """ Calculates the geodeticlength and thickness of a particle shape based on its area and diameter
    
    Calculation according to DIN ISO 9276-6: The geodetic lengths and thickness are approximated by an rectangle
    with the same area and perimeter:
    area = geodeticlength * thickness
    perimeter 2 * (geodeticlength + thickness)   
    
    Args:
        area (float): The area of the particle
        perimeter (float): The perimeter of the particle
        
    Returns:
        geodeticlength (float): The geodetic length of the particle
        thickness (float): The thickness of a particle
    
    """
    # helping variable: value under square root (see below)
    v = (perimeter/4)**2 - area
    
    # make sure value under squareroot is > 0
    v = max(v, 0)
    
    geodeticlength = perimeter/4 + np.sqrt(v)
    thickness = perimeter/2 - geodeticlength
    
    return geodeticlength, thickness
    

def calc_min_bounding_box(contour):    
    """ Calculates the minimal bounding box of an image based on a contour list.
    
    Args:
        contour (ndarray): List of contour points describing the particle in the form [[x0, y0], [x1, y1], ...].
        
    Returns:
        width_min_bb (float): The width of the minimal bounding box.
        height_min_bb (float): The height of the minimal bounding box.
            (height_min_bb >= width_min_bb, height_min_bb orthogonal to width_min_bb)
        center_bb (ndarray, float): Center of the minimal bounding box in the form [x,y]
        cornerpoints_min_bb (ndarray, float): A list of the four cornerpoints of the minimal bounding box.
            In form: [[x0, y0], [x1, y1], [x2, y2], [x3, y3]]
    """
    
    # find minimal bounding box
    min_bb_rect = cv2.minAreaRect(contour)
    
    # extract width, height and center
    (center_bb_y,center_bb_x), (height_min_bb, width_min_bb), alpha_bb = min_bb_rect
    center_min_bb = np.zeros((2,))
    center_min_bb[0], center_min_bb[1] = center_bb_x, center_bb_y
    
    # extract corner points
    box_cv2 = cv2.boxPoints(min_bb_rect)
    cornerpoints_min_bb = np.zeros((4,2))
    cornerpoints_min_bb[:,0], cornerpoints_min_bb[:,1] = box_cv2[:,0], box_cv2[:,1] #cv2 seems not to be consistent here
    
    return width_min_bb, height_min_bb, center_min_bb, cornerpoints_min_bb   

def calc_convex_perimeter(bw_convex):
    """ Calculates the perimeter of the convex hull based on an binary image of the convex hull.
    
    Args:
        bw_convex (ndarray, boolean): Binary image of the convex hull.
        
    Returns:
        convex_perimeter (float): Perimeter of the convex hull.
    """
    labels = measure.label(bw_convex)
    props = measure.regionprops_table(measure.label(bw_convex), properties=('label','perimeter'))
    convex_perimeter = props['perimeter'][0]
    return convex_perimeter

def calc_n_erosions_to_erase_binary_img(bw, pad_width=1, selem=None, return_list=False):
    """ Determines the number of erosions that are necessary to fully erase all true elements in a binary image.
    
    Args:
        bw (ndarray, boolean): Binary image.
        pad_width (int): Number of false pixels to pad around the image.
            (If outer pixels are not false, this effects the number of erosions (depneing on the neighborhood))
        selem (ndarray): The neighborhood expressed as a 2-D array of 1’s and 0’s.
            If None, use a cross-shaped structuring element (connectivity=1).
            (this parameter is directly passed to skimage.morphology.binary_erosion)
        return_list: If True, a list with the number of true pixels for each iteration is passed,
            if False, only the number of erosion (length of list) is returned.
            
    Returns:
        n_erosions (int) [if return_list == True]: Number of erosions to fully erase true pixels from image.
        n_true_pixels_list (list) [if return_list == False]: a list with the number of true pixels for each iteration.    
    """
    
    # Apply padding
    bw_eroded = np.pad(bw, pad_width=pad_width, mode='constant', constant_values=False)
    
    n_true_pixels_list = []
    n_true_pixels_list.append(np.count_nonzero(bw_eroded))

    while n_true_pixels_list[-1] > 0:
        bw_eroded = morphology.binary_erosion(bw_eroded)
        n_true_pixels_list.append(np.count_nonzero(bw_eroded))
        
    n_erosions = len(n_true_pixels_list)
        
    if return_list == True:
        return n_true_pixels_list
    else:
        return n_erosions
    
    
def calc_n_erosions_to_erase_binary_complement(bw_cropped, bw_convex_cropped, pad_width=1, selem=None, return_list=False):
    """ Number of erosions that are necessary to erase all pixels from the complement between the binary image of
        a particle and the binary image of its convex hull.
        
    Args:
        bw (ndarray, boolean): Binary image of particle.
        bw_convex_hull (ndarray, boolean): Binary image of the convex hull of the particle.
        pad_width (int): Number of false pixels to pad around the image.
            (If outer pixels are not false, this effects the number of erosions (depneing on the neighborhood))
        selem (ndarray): The neighborhood expressed as a 2-D array of 1’s and 0’s.
            If None, use a cross-shaped structuring element (connectivity=1).
            (this parameter is directly passed to skimage.morphology.binary_erosion)
        return_list: If True, a list with the number of true pixels for each iteration is passed,
            if False, only the number of erosion (length of list) is returned.
        
    Returns:
        n_erosions (int) [if return_list == True]: Number of erosions to fully erase true pixels from image.
        n_true_pixels_list (list) [if return_list == False]: a list with the number of true pixels for each iteration.    
    """
    assert bw_cropped.shape == bw_convex_cropped.shape, "Binary image and complement should have same shape, but have shape bw: {} and shape bw_complement: {}".format(bw.shape, bw_convex_hull.shape)
    
    bw_non_object = np.logical_not(bw_cropped)
    bw_complement = np.logical_and(bw_convex_cropped, bw_non_object)

    result = calc_n_erosions_to_erase_binary_img(bw_complement, pad_width=pad_width, selem=selem, return_list=return_list)
    return result


def first_nonzero(arr, axis, invalid_val=-1):
    mask = arr!=0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)


def calc_feret_diameter(bw):
    """ Calculates the Feret diameter of a particle (orthogonal to y-direction).

    Args:
        bw (ndarray, bool): Binary image.
        
    Returns:
        feret_diameter (int): Feret diameter.
        x_idx (ndarray, int): x-Coordinates of the two Feret calipers [x0, x1].
    """
    
    n_pixels_in_xdirection  = bw.shape[0]
    
    # "squeeze" bw in y-direction into one column
    n_true_pixels_in_ydirection = np.count_nonzero(bw, axis=1)
    
    # caliper from top    
    idx_first_nonzero_pixel = first_nonzero(n_true_pixels_in_ydirection, axis=0)
    if idx_first_nonzero_pixel == -1:
        # no element found
        return 0
    
    # caliper from bottom
    n_true_pixels_y_reversed = n_true_pixels_in_ydirection[::-1]
    idx_first_nonzero_pixel_from_bottom = first_nonzero(n_true_pixels_y_reversed, axis=0)
    idx_last_nonzero_pixel = n_pixels_in_xdirection - idx_first_nonzero_pixel_from_bottom
    
    # feret diameter
    feret_diameter = idx_last_nonzero_pixel - idx_first_nonzero_pixel
    x_idx = np.zeros((2,))
    x_idx[0] = idx_first_nonzero_pixel
    x_idx[1] = idx_last_nonzero_pixel
    return feret_diameter, x_idx


def calc_martin_diameter(bw, area_share_from_bottom=0.5, return_index=False):
    """ Calculates the Martin diameter of a particle (in y-direction).
    
    Note: The original Martin diameter is is measured at the x-position, where the particle is split into 50% / 50% area share.
        However, we can also calculate a diameter at a flexible area share. This is given by area_share_from_bottom.
        area_share_from_bottom = 1 - area_share_from_top

    Args:
        bw (ndarray, bool): Binary image.
        area_share_from_bottom (float, optional): area share, where the Martin diameter is measured.
            Default is 0.5, the original definition of the Martin diameter.
        
    Returns:
        martin_diameter (int): Martin diameter.
        idx (ndarray, int): Indexes of the start (0) and endpoint (1) of the Martin diameter in form [[x0,y0], [x1,y1]].
    """
    if area_share_from_bottom == 0 or area_share_from_bottom == 1:
        return 0
    elif area_share_from_bottom < 0 or area_share_from_bottom > 1:
        print("Invalid area share.")
        return None    
    
    area_share_from_top = 1 - area_share_from_bottom
    
    # calculate the number of True pixels in ydirection
    n_true_pixels_in_ydirection = np.count_nonzero(bw, axis=1)
    # get the index where the area_share is reached (from top)
    cum_area = np.cumsum(n_true_pixels_in_ydirection)
    area = cum_area[-1]
    if area == 0:
        print("No object found.")
        return 0
    greater_than_area_share = cum_area > (area_share_from_top * area)
    x_split = greater_than_area_share.argmax()

    # extract row of this position
    row_split = bw[x_split,:]
    
    # calculate the martin diameter of this dimension
    n_pixels = row_split.shape[0]

    first_y_idx = row_split.argmax()
    row_split_reversed = row_split[::-1]
    last_y_idx = n_pixels - row_split_reversed.argmax()

    martin_diameter = last_y_idx - first_y_idx
        
    idx = np.zeros((2,2))
    idx[0,0], idx[0,1] = x_split, first_y_idx
    idx[1,0], idx[1,1] = x_split, last_y_idx
    
    return martin_diameter, idx


def calc_longest_chord(bw):
    """ Calculates the maximum chords of a particle shape in y-direction.
    
    Args:
        bw (ndarray, bool): Binary image.
        
    Returns:
        max_chord (int): Max chord.
        line_max_chord (ndarray, int): Start (a) and endpoint (b) of max chord in form [[xa,ya],[xb,yb]]
    """
    
    all_chords, all_edgepoints = calc_chords(bw)
    max_chord, line_max_chord = determine_max_chord(all_chords, all_edgepoints)
    
    return max_chord, line_max_chord


def determine_max_chord(all_chords, all_edgepoints):
    """ Determines the maximum of all chords.
    
    Args:
        all_chords (ndarray, int): Array of several chords in form [c0, c1, ...]
        all_edgepoints (ndarray, int): List of edgepoints referrring to the chords
            in form [[x0a,y0a],[x0b,y0b],[x1a,y1a],[x1b,y1b],...], where a is the start
            and b is the endpoint of the corresponding edge.
            
    Returns:
        max_chord (int): Max chord.
        line_max_chord (ndarray, int): Start (a) and endpoint (b) of max chord in form [[xa,ya],[xb,yb]]
    """
    max_chord = all_chords.max()
    idx = all_chords.argmax()
    idx_points = 2 * idx
    
    line_max_chord = all_edgepoints[idx_points:idx_points+2, :]
    
    return max_chord, line_max_chord



def calc_chords(bw):
    """ Calculates all chords of a particle shape in y-direction.
    
    Args:
        bw (ndarray, bool): Binary image.
        
    Returns:
        chords (ndarray, int): An array with length of all found chords in y-direction.
    """
    # Initialize array to store found chords
    # Note: There may be none or several (> 1) choords per row (this is why we use the append command)
    all_chords = np.array([], dtype='int64')
    all_edgepoints = np.zeros((0,2), dtype='int64')

    # Find points where values changes in y-direction (horizontal)
    # So from False to True or True to False (i.e. two neighbor pixels in y-direction have different values)
    # Background pixels are False and Object Pixels are True. To make sure we find also the first and last
    # changing points, we add a false padding to the first and last column

    # padding
    bw_pad = np.pad(bw, pad_width=1, mode='constant', constant_values=False)
    
    bw_rep = np.repeat(bw_pad, repeats=1, axis=1)
    
    # find changing points in y-direction
    idx_points_with_changes = np.array(np.where(bw_rep[:,:-1] != bw_rep[:,1:]))

    # Reshape the results
    x_idx = idx_points_with_changes[0,:]
    y_idx = idx_points_with_changes[1,:]
    points = np.column_stack((x_idx, y_idx))

    # we are only interested in the rows, where changes accur => get indexes of these rows
    unique_x_positions = np.unique(x_idx)
    
    
    for u in unique_x_positions:
        # extract points of this row
        points_this_row = points[np.where(points[:,0]==u)]

        # to calculate the chords we now need to determine the distance between pairs of changepoints in y-direction
        y_coords_this_row = points_this_row[:,1]
        
        # since the leftest row is always filled with False pixels (see padding above), we know that the first point
        # will be a changing point from False to True, the next point then will be a changing point from True to
        # False and so on. More generally: If we sort the y-coordinates of all changing-points in ascending order and or
        # indexes start at zero, then: All changing points with even indexes (0, 2, ...) are starting points
        # and all ending points with uneven indexes (1, 3, ...) are ending points

        # as promised: sort y_coordinates in ascending order
        y_coords_this_row_sorted = np.sort(y_coords_this_row)        
         
        # starting points have even indexes
        starting_points_this_row = y_coords_this_row_sorted[0::2]

        # ending points have uneven indexes
        end_points_this_row = y_coords_this_row_sorted[1::2]

        # The choords are the distance between the corresponding starting and ending points
        chords_this_row = end_points_this_row - starting_points_this_row
        # add found choords to result list
        all_chords = np.append(all_chords, chords_this_row)
        
        # undo padding for idx
        points_this_row = points_this_row - 1
        
        all_edgepoints = np.concatenate((all_edgepoints, points_this_row))
        
        # undo padding shift:
       # all_edgepoints = all_edgepoints - 1
        
    return all_chords, all_edgepoints


def calc_nassenstein_diameter(bw):
    
    """ Calculates the Nassenstein diameter of a particle shape.
    
    Note: There might be several touching points in the lowest row.
        In this implementation we will evaluate the Nassenstein Durchmesser at the middle
        of the continuous first contact surface from left.
    
    Args:
        bw (ndarray, bool): Binary image.
        
    Returns:
        nassenstein_diameter (int): Nassenstein diameter.
        idx (ndarray, int): Indexes of the start (0) and endpoint (1) of the Nassenstein diameter in form [[x0,y0], [x1,y1]].
    """
    # padding
    bw_pad = np.pad(bw, pad_width=1, mode='constant', constant_values=False)
    
    # find lowest row
    n_pixels_xdirection = bw_pad.shape[0]
    n_true_pixels_in_ydirection = np.count_nonzero(bw_pad, axis=1)
    n_true_pixels_in_ydirection_from_bottom = n_true_pixels_in_ydirection[::-1]
    idx_lowest_row = n_pixels_xdirection - first_nonzero(n_true_pixels_in_ydirection_from_bottom,axis=0)
    idx_lowest_row = idx_lowest_row - 1
    lowest_row = bw_pad[idx_lowest_row, :]
    # obtain first touching surface by finding the first two changing points:
    changing_points_row = np.where(lowest_row[:-1] != lowest_row[1:])[0]
    changing_points_row = np.sort(changing_points_row)

    start_idx_first_contact = changing_points_row[0]
    end_idx_first_contact = changing_points_row[1]

    evaluation_idx = int((start_idx_first_contact + end_idx_first_contact)/2)

    # extract the column at evaluation idx:
    nassenstein_column = bw_pad[:,evaluation_idx]

    # we measure starting from bottom
    nassenstein_column_from_bottom = nassenstein_column[::-1]

    # again we consider the changing points to determine the Nassenstein diameter
    changing_idx_nassenstein_column = np.where(nassenstein_column_from_bottom[:-1] != nassenstein_column_from_bottom[1:])[0]
    changing_idx_nassenstein_column = np.sort(changing_idx_nassenstein_column)
    

    # since we started counting from bottom, we have to transform the indexes to the coordinate system from top
    measurement_point_bottom = n_pixels_xdirection - changing_idx_nassenstein_column[0]
    measurement_point_top = n_pixels_xdirection - changing_idx_nassenstein_column[1]

    nassenstein_diameter = measurement_point_bottom - measurement_point_top
    
    idx = np.zeros((2,2), dtype=('int64'))
    idx[0,0], idx[0,1] = measurement_point_top, evaluation_idx
    idx[1,0], idx[1,1] = measurement_point_bottom, evaluation_idx
    
    # undo padding shift:
    idx = idx - 1

    return nassenstein_diameter, idx


def calc_distribution_parameters(statistical_length):
    
    max_value = statistical_length.max()
    min_value = statistical_length.min()
    median_value = np.median(statistical_length)
    mean_value = np.median(statistical_length)
    mode, counts = stats.mode(statistical_length)
    std = np.std(statistical_length)
    
    return max_value, min_value, median_value, mean_value, mode, std
    

def calc_statistical_length_distribution(bw, daplha):
    """ Calculates the statistical length (Feret-, Martin-, Nassenstein-diameter, chords and max chord) for a
        binary image in daplha degree steps.
        
        Args:
        bw (ndarray, bool): Binary image.
        dalpha (int/float): Rotation stepsize in degree (0 - 180°).
        
        Returns:
            feret_diameters (ndarray, int): Array of Feret diameters for each rotation.
            martin_diameters (ndarray, int): Array of Martin diameters for each rotation.
            nassenstein_diameters (ndarray, int): Array of Nassenstein diameters for each rotation.
            max_chords (ndarray, int): Array of maximum chord for each rotation.
            all_chords (ndarray, int): Array of all chords for each rotation.
            measured_angles (ndarray, float): Array of the rotated angles determinated by daplha.
    """

    angles = np.arange(0,180,daplha)
    # introduce empty arrays
    feret_diameters = np.zeros(angles.shape, dtype='int64')
    martin_diameters = np.zeros(angles.shape, dtype='int64')
    nassenstein_diameters = np.zeros(angles.shape, dtype='int64')
    max_chords = np.zeros(angles.shape, dtype='int64')
    all_chords = np.array([])

    # iterate over all angles
    for i, angle in enumerate(angles):
        # rotate image
        bw_rotated = transform.rotate(bw, angle, resize=True)
        # important: skimage.transform.rotate returns a greyscale image with values between 0 and 1 
        # => use simple thresholding to transform back to binary image
        bw_rotated = bw_rotated > 0.5

        # Feret diameter
        feret_diameters[i], _ = calc_feret_diameter(bw_rotated)

        # Martin diameter
        martin_diameters[i], _ = calc_martin_diameter(bw_rotated)

        # Nassenstein diameter
        nassenstein_diameters[i], _ = calc_nassenstein_diameter(bw_rotated)

        # Chords
        chords, _ = calc_chords(bw_rotated)
        max_chords[i] = chords.max()
        all_chords = np.append(all_chords, chords)
        
        measured_angles = angles
    return feret_diameters, martin_diameters, nassenstein_diameters, max_chords, all_chords, measured_angles 



def max_dimensions(max_chords, angles):
    """ Calculates the max dimensions of a particle shape.
    
    x_max is the overall max chord of the particle in all possible orientations. y_max is the longest chord orthogonal to
    y_max.
    
    Args:
        max_chords (ndarray, int/float): Array of max_chords of a particle at different angles in shape [c0, c1, c2, ...].
        angles (ndarray, int/float): The respective angles to the max_chords in ascending order.
    
    Returns:
        x_max (int/float): Larger max dimension of the particle (definition see above).
        y_max (int/float): Smaller max dimension of the particle (definition see above).
    """
    
    assert max_chords.shape == angles.shape, "max_chords and angles should have the same shape."
    
    
    x_max = max_chords.max()
    idx_x_max = max_chords.argmax()
    angle_x_max = angles[idx_x_max]
    
    angle_y_max = (angle_x_max + 90) % 180
    idx_y_max = (np.abs(angles - angle_y_max)).argmin()
    y_max = max_chords[idx_y_max]
    
    return x_max, y_max 


def approximate_fractal_dimension_by_slope(measurements, step_sizes):
    """ Approximates fractal dimension form measurements and step_sizes by slope the log/log Richardson plot.
    
            Slope is approximate by linear regression of the curve in the log/log Richardson plot.
            
        Args:
            measurements (ndarray, float): Array of measured length of corresponding step_sizes.
            step_sizes (ndarray, float/int): Correspong step_sizes. Must have same shape as measurements
            
        
        Returns:
            DF (float): Approximated fractal dimension.
    """
    
    assert measurements.shape == step_sizes.shape, "Measurements and step size must have same shape."

    # remove entries where measurements is zero
    measurements_org = measurements.copy()
    measurements = measurements[measurements_org != 0]
    step_sizes = step_sizes[measurements_org != 0]
    
    measurements_log2 = np.log2(measurements)
    step_sizes_log2 = np.log2(step_sizes)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(step_sizes_log2, measurements_log2)
    
    DF = 1 - slope
    
    return DF


def uniformly_structured_gait(contour, step_sizes):
    """ Performs a uniformly_structured_gait of different stepsizes and returns the walked perimeter.
    
    Args:
        contour (ndarray, int/float): Array of contourpoints describing the shape contour
            in form [[x0, y0], [x1, y1], [x2, y2]]
        step_sizes (ndarray, float): Array of the different step sizes to take.
            
    Returns:
        perimeters(ndarray, float): Array of walked perimeters in shape [p0, p1, p2, ...].
    """

    perimeters = np.zeros(step_sizes.shape)

    all_idx_list = []

    for j, step_size in enumerate(step_sizes):
        idx_list = [0]
        distance_list = []

        start_point = contour[0]
        i = 1
        while i < len(contour):
            end_point = contour[i]
            dist = spatial.distance.euclidean(start_point, end_point)
            if dist >= step_size:
                idx_list.append(i)
                distance_list.append(dist)
                start_point = contour[i]
                i = i + 1

            i = i + 1

        all_idx_list.append(idx_list)

        perimeters[j] = sum(distance_list)
    
    return perimeters


def calc_fractal_dimension_perimeter_method(contour, max_feret, n_stepsizes=10, step_size_min=2):
    """ Approximates the fractal dimension based on the perimeter method.
            Note: Stepsizes range from step_size_min to step_size_max in a log2 grid.
            step_size_max is set to 0.3 * max_feret according to DIN ISO 9276-6.
            
        Args:
            contour (ndarray, int/float): Array of contourpoints describing the shape contour
                in form [[x0, y0], [x1, y1], [x2, y2]]
            max_feret (float): Max feret diameter of the particle shape for norming the perimeters.
            n_stepsizes (int): Number of different stepsizes to take. Default is 10.
            step_size_min (float): Minimum stepsize to walk. Default is 1. (Definition of max stepsize see above)
            
        
        Returns:
            DF (float): Approximated fractal dimension.
    """


    step_size_max = 0.3 * max_feret
    step_sizes = np.logspace(np.log2(step_size_min), np.log2(step_size_max), num=n_stepsizes, endpoint=True, base=2)

    perimeters = uniformly_structured_gait(contour, step_sizes)
    
    # Normalize by maximum feret diameter
    perimeters_normed = perimeters/max_feret
    
    # Determine fractal dimension
    DF = approximate_fractal_dimension_by_slope(perimeters_normed, step_sizes)

    return DF



def calc_fractal_dimension_boxcounting_method(bw, pad_width=1):
    """ Approxiamte the fractal dimension of a binary image by the boxcounting method.
    
    Args:
        bw (ndarray, bool): Binary image.
        pad_width (int, optional): Width of applied zero-padding around the image. Default is 1.
        
    Returns:
        DF (float): The approximated fractal dimension.
    
    """
    bw = np.pad(bw, pad_width=pad_width, mode='constant', constant_values=False)
    number_of_boxes, box_sizes = box_counting(bw)
    DF = approximate_fractal_dimension_by_slope(number_of_boxes, box_sizes)
    
    return DF
    

def box_counting(bw, min_box_size=2):
    """ Counts the number of non-zero and non-full boxes of and binary image at different box sizes.
    
        Args:
            bw (ndarray, bool): Binary image.
            min_box_size (int): Minimum investigated box size. Default is 2. Must be representable as 2**n, where 
                n is a natural number.
            
        Returns:
            number_of_boxes (ndarray, int): Array if the numbers of found non-zero and non-full boxes for the corresponding box size.
            box_sizes (ndarray, int): Array of corresponding box sizes.
    """

    # make bw shape (2**n, 2**n)
    bw_pad = zeropad_bw_to_shape_of_power_two(bw)

    # create box sizes
    bw_size = bw_pad.shape[0]
    max_box_size = bw_size/2
    exp_max_box, exp_min_box = int(np.log2(max_box_size)), int(np.log2(min_box_size))
    n_steps = exp_max_box - exp_min_box + 1
    box_sizes = np.logspace(exp_min_box, exp_max_box, num=n_steps, base=2, dtype='int', endpoint=True)
    
    # initialize array for solutions
    number_of_boxes = np.zeros(box_sizes.shape, dtype='int64')

    # determine number of boxes for different box sizes
    for i, box_size in enumerate(box_sizes):
        number_of_boxes[i] = boxcount_single_boxsize(bw_pad, box_size)
        
    return number_of_boxes, box_sizes


def zeropad_bw_to_shape_of_power_two(bw):
    """ Transform a binary image into shape (2**n, 2**n) (where n is a natural number) by placing it on
        a black background.
        
        Args:
            bw (ndarray, bool): Binary image.
            
        Returns:
            bw_pad (ndarray, bool): Padded binary image of shape (2**n, 2**n), where n is a natural number.
    """

    max_bw_shape = max(bw.shape)

    # determine shape of the padded image
    exponent_bw_shape = int(np.ceil(np.log2(max_bw_shape)))
    padded_bw_size = 2**exponent_bw_shape
    padded_bw_shape = (padded_bw_size, padded_bw_size)

    # initialize the padded image
    bw_pad = np.zeros(padded_bw_shape, dtype='bool')

    # determine shift, s.t. bw is inserted in the center of bw_pad
    shift_x, shift_y = int((padded_bw_size - bw.shape[0])/2), int((padded_bw_size - bw.shape[1])/2)

    # insert bw
    bw_pad[shift_x:bw.shape[0]+shift_x, shift_y:bw.shape[1]+shift_y] = bw
    
    return bw_pad


def boxcount_single_boxsize(bw, box_size):
    """ Calculates the number of boxes of shape (box_size, box_size) that are non-empty and non-full in an image.

    Args:
        bw (ndarray, bool): Binary image.
        box_size (int): Size of the boxes. Boxes a value of 2**n (n in natural numbers)

    Returns:
        n_boxes (int): The number of found boxes.
    """


    # From [https://github.com/rougier/numpy-100 (#87)] cited from [https://gist.github.com/viveksck/1110dfca01e4ec2c608515f0d5a5b1d1]
    # create boxes of shape (box_size,box_size) and store results in a matrix
    # an efficient implementation of this procedere is:
    S = np.add.reduceat(np.add.reduceat(bw, np.arange(0, bw.shape[0], box_size), axis=0),
                        np.arange(0, bw.shape[1], box_size), axis=1)

    # count non-empty (0) and non-full boxes (box_sizes*box_size)
    n_boxes = len(np.where((S > 0) & (S < box_size*box_size))[0])
    return n_boxes


