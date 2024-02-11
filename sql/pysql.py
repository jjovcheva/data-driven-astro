import psycopg2
import numpy as np

def select_all(tablename):
    '''Query table and return all rows.'''
    
    # Establish the connection.
    conn = psycopg2.connect(dbname='db', user='mysqlusername')

    # Create cursor object to interface with the database.
    cursor = conn.cursor()
    
    cursor.execute('select * from %s;' %tablename)
    return cursor.fetchall()

def column_stats(tablename, columnname):
    '''Calculate mean and median of a selected column.'''
    
    # Establish the connection.
    conn = psycopg2.connect(dbname='db', user='mysqlusername')

    # Create cursor object to interface with the database.
    cursor = conn.cursor()
    
    cursor.execute('select %s from %s;' %(columnname, tablename))
    records = cursor.fetchall()
    array = np.array(records)
    
    return np.mean(array), np.median(array)

def query_rad(filename):
    '''
    Replicate a basic SQL query that returns the kepler ID
    and radius of stars larger than the Sun, in ascending radius.
    '''
    # Read in data.
    data = np.loadtxt(filename, delimiter=',', usecols=[0, 2])
    
    # Examine radii and filter to keep only those larger than the Sun.
    data = data[data[:, 1] > 1, :]
    
    # Sort in ascending order of radius.
    sorting = np.argsort(data[:, 1])
    data = data[sorting, :]
    return data

def query_ratio(s_file, p_file):
    '''
    Replicate a basic SQL query that finds the ratio between the radii
    of stars larger than our Sun and their corresponding planets, and
    returns the ratios in ascending order.
    '''
    stars = np.loadtxt(s_file, delimiter=',', usecols=[0, 2])
    planets = np.loadtxt(p_file, delimiter=',', usecols=[0, 5])
    
    # Examine radii and filter to keep only those larger than the Sun.
    s_data = stars[stars[:, 1] > 1, :]
    
    # Create a list to store all radius ratios.
    r_ratios = []

    # Iterate over each filtered star.
    for star in s_data:
        # Find corresponding planets for each star.
        p_matching = planets[planets[:, 0] == star[0]]
        
        for planet in p_matching:
            # Calculate radius ratio and append to the list.
            r_ratio = planet[1] / star[1]
            r_ratios.append(r_ratio)

    r_ratios = np.array(r_ratios)
    r_ratios = r_ratios.reshape(-1, 1)
    
    # Sort in ascending order of radius ratio.
    sorting = np.argsort(r_ratios[:, 0])
    return r_ratios[sorting]
