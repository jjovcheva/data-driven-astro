import numpy as np
from astropy.io import fits
import time
from astropy.coordinates import SkyCoord
from astropy import units as u

def calc_stats(file_path):
    """
    Calculate the mean and median of a 1D array stored in a CSV file.

    Parameters
    ----------
    file_path : str
        Path to the CSV file containing the 1D array.

    Returns
    -------
    tuple
        A tuple containing the mean and median of the array.
    """
    data = np.loadtxt(file_path, delimiter=',')
    return np.round(np.mean(data), 1), np.round(np.median(data), 1)

def mean_datasets(filenames):
    """
    Calculate the mean of a set of signals stored in CSV files.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing the signals.

    Returns
    -------
    numpy.ndarray
        The mean of the signals.
    """
    n = len(filenames)
    data = np.loadtxt(filenames[0], delimiter=',')
    for filename in filenames[1:]:
        data += np.loadtxt(filename, delimiter=',')
    
    mean = data / n
    return np.round(mean, 1)

def brightest(filename):
    """
    Open FITS file and find the coordinates of the brightest pixel.

    Parameters
    ----------
    filename : str
        Path to the FITS file.

    Returns
    -------
    tuple
        Coordinates (row, column) of the brightest pixel.
    """
    hdulist = fits.open(filename)
    data = hdulist[0].data
    hdulist.close()
    return np.unravel_index(data.argmax(), data.shape)

def mean_fits(filenames):
    """
    Calculate the mean of each pixel over a set of FITS files.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing FITS files.

    Returns
    -------
    mean: numpy.ndarray
        The mean of each pixel.
    """
    n = len(filenames)
    hdulist = fits.open(filenames[0])
    data = hdulist[0].data
    hdulist.close()

    for filename in filenames[1:]:
        hdulist = fits.open(filename)
        data += hdulist[0].data
        hdulist.close()

    mean = data / n
    return mean

def time_stat(func, size, ntrials):
    """
    Measure the time required for different functions to process a random array.

    Parameters
    ----------
    func : function
        The function to be timed.
    size : int
        Size of the random array to be processed.
    ntrials : int
        Number of trials for timing.

    Returns
    -------
    float
        The average time taken for the function over multiple trials.
    """
    t_tot = 0
    for _ in range(ntrials):
        data = np.random.rand(size)
        start = time.perf_counter()
        _ = func(data)
        t_tot += time.perf_counter() - start

    return t_tot / ntrials

def median_fits(filenames):
    """
    Calculate the median of each pixel over a set of FITS files.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing FITS files.

    Returns
    -------
    numpy.ndarray
        The median of each pixel.
    float
        Time taken for the operation.
    float
        Memory usage in kilobytes.
    """
    start = time.time()
    data = [fits.open(file)[0].data for file in filenames]
    median = np.median(np.dstack(data), axis=2)
    memory = np.dstack(data).nbytes / 1024
    end = time.time() - start
    return median, end, memory

def median_bins(values, nbins):
    """
    Split list of values into bins to calculate the median.

    Parameters
    ----------
    values : list or numpy.ndarray
        List or array of floats.
    nbins : int
        Number of bins to split the values.

    Returns
    -------
    tuple
        Tuple containing mean, standard deviation, bin index of values < (mean - std),
        histogram of values within (mean - std) to (mean + std), and bin width.
    """
    mean = np.mean(values)
    std_dev = np.std(values)
    
    # Bin holding values less than (mean - std).
    min_bin = 0
    bins = np.zeros(nbins)
    
    minval = mean - std_dev
    maxval = mean + std_dev
    width = 2 * std_dev / nbins
    
    for val in values:
      if val < minval:
        min_bin += 1
      elif val < maxval:
        bin = int((val - minval)/width)
        bins[bin] += 1
                  
    return mean, std_dev, min_bin, bins, width
  
def median_approx(values, nbins):
    """
    Use binapprox algorithm to calculate the median of a list.

    Parameters
    ----------
    values : list or numpy.ndarray
        List or array of floats.
    nbins : int
        Number of bins to split the values.

    Returns
    -------
    float
        Median of the values.
    """
    mean, std_dev, min_bin, bins, width = median_bins(values, nbins)
    
    N = len(values)
    mid_val = np.multiply(N+1, 0.5)

    count = min_bin 
    for b, bincount in enumerate(bins):
      count += bincount
      if count >= mid_val:
        break

    median = mean - std_dev + np.multiply(width, b+0.5)

    return median

def update_stats(filenames):
    """
    Use Welford's algorithm to update the mean and standard deviation of a set of FITS files.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing FITS files.

    Returns
    -------
    tuple
        Tuple containing the updated mean and standard deviation.
    """
    count, mean, sigma = 0, 0, 0
    for filename in filenames:
      hdu_list = fits.open(filename)
      data = hdu_list[0].data
      
      count += 1
      delta = data - mean
      mean += delta / count
      sigma += delta * (data-mean)
      hdu_list.close()
      
    if count < 2:
      return mean, None
    else:
      mean, var, sample_var = mean, np.divide(sigma, count), np.divide(sigma, count-1)
      
    std_dev = np.sqrt(sample_var)
    return mean, std_dev               

def median_bins_fits(filenames, nbins):
    """
    Split array elements from FITS file into bins.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing FITS files.
    nbins : int
        Desired number of bins to split files into.

    Returns
    -------
    tuple
        Tuple containing mean, standard deviation, bin index of values < (mean - std),
        histogram of values within (mean - std) to (mean + std), and bin width.
    """
    mean, std = update_stats(filenames)
    dim = mean.shape # Dimension of the FITS file arrays

    # Initialise bins
    left_bin = np.zeros(dim)
    bins = np.zeros((dim[0], dim[1], nbins))
    bin_width = 2 * std / nbins

    # Loop over all FITS files
    for filename in filenames:
        hdulist = fits.open(filename)
        data = hdulist[0].data

        # Loop over every point in the 2D array
        for i in range(dim[0]):
          for j in range(dim[1]):
            value = data[i, j]
            m = mean[i, j]
            s = std[i, j]
            
            minval = m - s
            maxval = m + s

            if value < minval:
              left_bin[i, j] += 1
                  
            elif value >= minval and value < maxval:
              bin = int((value - minval)/bin_width[i, j])
              bins[i, j, bin] += 1

    return mean, std, left_bin, bins

def median_approx_fits(filenames, nbins):
    """
    Calculate median of all elements across file stack.

    Parameters
    ----------
    filenames : list of str
        List of file paths containing FITS files.
    nbins : int
        Desired number of bins to split files into.

    Returns
    -------
    numpy.ndarray
        Median of each element across the file stack.
    """
    mean, std, left_bin, bins = median_bins_fits(filenames, nbins)
      
    dim = mean.shape # Dimension of the FITS file arrays
      
    N = len(filenames)
    mid_val = np.multiply(N+1, 0.5)
    bin_width = 2 * std / nbins
  
    median = np.zeros(dim)   
    
    for i in range(dim[0]):
      for j in range(dim[1]):    
        count = left_bin[i, j]
        for b, bincount in enumerate(bins[i, j]):
          count += bincount
          if count >= mid_val:
            break
        median[i, j] = mean[i, j] - std[i, j] + bin_width[i, j]*(b + 0.5)
        
    return median
      
def hms2dec(h, m, s):
    """
    Convert HMS coordinates to decimal degrees.

    Parameters
    ----------
    h : float
        Hours.
    m : float
        Minutes.
    s : float
        Seconds.

    Returns
    -------
    float
        Coordinate in decimal degrees.
    """
    m = np.divide(m, 60)
    s = np.divide(s, 3600)
    data = np.array([h, m, s])
    
    deg = float(np.multiply(15, np.sum(data)))
    return deg

def dms2dec(d, m, s):
    """
    Convert DMS coordinates to decimal degrees.

    Parameters
    ----------
    d : float
        Degrees.
    m : float
        Minutes.
    s : float
        Seconds.

    Returns
    -------
    float
        Coordinate in decimal degrees.
    """
    m = np.divide(m, 60)
    s = np.divide(s, 3600)
    
    if d < 0:
      data = np.sum(np.array([d*-1, m, s]))
    else:
      data = np.sum(np.array([d, m, s]))
      
    deg = float(np.where(d < 0, -1*data, data))
    
    return deg

def angular_dist(a_1, d_1, a_2, d_2, unit):
    """
    Calculate distance between 2 points on the celestial sphere.

    Parameters
    ----------
    a_1, d_1 : float
        RA and Dec of point 1.
    a_2, d_2 : float
        RA and Dec of point 2.
    unit : str
        Unit of the angular coordinates ('deg' or 'rad').

    Returns
    -------
    float
        Distance between 2 points.
    """
    if unit == 'deg':
      a_1, d_1, a_2, d_2 = np.radians(a_1), np.radians(d_1), np.radians(a_2), np.radians(d_2)
      a = np.sin(np.multiply(np.abs(d_1-d_2), 0.5)) ** 2
      b = np.cos(d_1) * np.cos(d_2) * np.sin(np.multiply(np.abs(a_1-a_2), 0.5)) ** 2
      
      d = 2 * np.arcsin(np.sqrt(a+b))
      d = np.degrees(d)
      
    elif unit == 'rad':
      a = np.sin(np.multiply(np.abs(d_1-d_2), 0.5)) ** 2
      b = np.cos(d_1) * np.cos(d_2) * np.sin(np.multiply(np.abs(a_1-a_2), 0.5)) ** 2
      
      d = 2 * np.arcsin(np.sqrt(a+b))
      
    return d

def import_bss(filename):
    """
    Read in bright radio source catalogue.
    Column 1 is the ID number.
    Columns 2-4 are the RA in HMS notation.
    Columns 5-7 are the Dec in DMS notation.

    Parameters
    ----------
    filename : str
        Path to the bright radio source catalogue file.

    Returns
    -------
    list
        List of tuples containing ID, RA in degrees, and Dec in degrees.
    """
  
    res = []
    cat = np.loadtxt(filename, usecols=range(1,7))
    for i, row in enumerate(cat, 1):
      res.append((i, hms2dec(row[0], row[1], row[2]), dms2dec(row[3], row[4], row[5])))
    
    return res

def import_super(filename):
    """
    Read in SuperCOSMOS all-sky catalogue (CSV format).
    Column 1 is the RA in decimal degrees.
    Column 2 is the Dec in decimal degrees.
    Column 3 is other data, e.g. magnitude.

    Parameters
    ----------
    filename : str
        Path to the SuperCOSMOS catalogue file.

    Returns
    -------
    list
        List of tuples containing ID, RA in degrees, and Dec in degrees.
    """
    res = []
    cat = np.loadtxt(filename, delimiter=',', skiprows=1, usecols=[0, 1])
    for i, row in enumerate(cat, 1):
      res.append((i, row[0], row[1]))
      
    return res

def find_closest(cat, ra, dec, unit):
    """
    Find the closest point in a catalog to a given RA and Dec.

    Parameters
    ----------
    cat : list
        Catalog containing tuples (ID, RA, Dec).
    ra, dec : float
        Target RA and Dec.
    unit : str
        Unit of the angular coordinates ('deg' or 'rad').

    Returns
    -------
    tuple
        Tuple containing the ID and angular distance of the closest point.
    """
    d_min = np.inf
    id_min = None
    
    for idx, ras, decs in cat:
      d = angular_dist(ras, decs, ra, dec, unit)
      if d < d_min:
        id_min = idx
        d_min = d
    
    return id_min, d_min
  
def crossmatch_simple(cat1, cat2, d_max, unit):
    """
    Perform a simple crossmatch between two catalogs.

    Parameters
    ----------
    cat1 : numpy.ndarray
        The first catalog.
    cat2 : numpy.ndarray
        The second catalog.
    d_max : float
        The maximum distance for a match.
    unit : str
        Unit of the angular coordinates ('deg' or 'rad').

    Returns
    -------
    list
        List of matches (source index, reference index, angular distance).
    list
        List of indices with no match.
    float
        Time taken for the crossmatch operation.
    """
    start = time.perf_counter()
    match = []
    no_match = []
    
    cat1 = np.radians(cat1)
    cat2 = np.radians(cat2)
    d_max = np.radians(d_max)
    
    cat2_sorting = np.argsort(cat2[:,1])
    cat2_sorted = cat2[cat2_sorting]
    dec2_sorted = cat2_sorted[:,1]
    
    ra2, dec2 = cat2[:,0], cat2[:,1]
    
    for id1, (ra1, dec1) in enumerate(cat1):
      d_min = np.inf
      id2_min = None 
      
      min_dec = dec1 - d_max
      max_dec = dec1 + d_max
      
      start_idx = np.searchsorted(dec2_sorted, min_dec, side='left')
      end_idx = np.searchsorted(dec2_sorted, max_dec, side='right')

      for id2, (ra2, dec2) in enumerate(cat2_sorted[start_idx:end_idx+1], start_idx):  
        d = angular_dist(ra1, dec1, ra2, dec2, unit)
        if d < d_min:
          id2_min = cat2_sorting[id2]
          d_min = d
      
      if d_min > d_max:
        no_match.append(id1)
      else:
        match.append((id1, id2_min, np.degrees(d_min)))
          
    time_taken = time.perf_counter() - start
    return match, no_match, time_taken  

def crossmatch(cat1, cat2, d_max):
    """
    Perform a crossmatch between two catalogs.

    Parameters
    ----------
    cat1 : numpy.ndarray
        The first catalog.
    cat2 : numpy.ndarray
        The second catalog.
    d_max : float
        The maximum distance for a match.

    Returns
    -------
    list
        List of matches (source index, reference index, angular distance).
    list
        List of indices with no match.
    float
        Time taken for the crossmatch operation.
    """
    start = time.perf_counter()
    match = []
    no_match = []
    
    sky_cat1 = SkyCoord(cat1*u.degree, frame='icrs')
    sky_cat2 = SkyCoord(cat2*u.degree, frame='icrs')
        
    closest_ids, closest_dists, _ = sky_cat1.match_to_catalog_sky(sky_cat2)
    for id1, (id2_min, d) in enumerate(zip(closest_ids, closest_dists)):
      d_min = d.value
      if d_min > d_max:
        no_match.append(id1)
      else:
        match.append([id1, id2_min, d_min])
          
    time_taken = time.perf_counter() - start
    return match, no_match, time_taken

