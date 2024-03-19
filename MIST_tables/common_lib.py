'Common functions'

import numpy as np
import os

def MIST_masses():
    """

    Creates the array of all de star masses of the MIST tables.

    Args:
        None.

    Usage:
        >> MIST_masses = MIST_masses()
    
    """
    mass01_04 = np.arange(10, 40, 2) / 100 #0.02 steps
    mass04_09 = np.arange(40, 90, 5) / 100 #0.05 steps
    mass09_11 = np.arange(90, 110, 2) / 100 #0.02 steps
    mass11_13 = np.arange(110, 135, 5) / 100 #0.05 steps
    mass = np.concatenate((mass01_04, mass04_09, mass09_11, mass11_13))
    return mass

def find_nearest(array, value):
    """

    Find nearest neighbour in array.

    Args:
        None.

    Usage:
        >> nearest_neighbour = find_nearest(array, value)
    
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_2_nearest(array, value):
    """

    This function finds the 2 array values nearest to the input value.
    Returns the 2 nearest values from the array.

    Args:
        array: array with some values from which you want to find the 2 nearest to the input value
        value: your input value

    Usage:
        >> value1, value2 = find_2_nearest(array, value)
    
    """

    array = np.asarray(array)
    dif = np.abs(array - value)
    sort = np.argsort(dif)
    idx = sort[0] #Nearest neighbour
    idx2 = sort[1] #2nd nearest neighbour
    return array[idx], array[idx2]

def interpolation(a, b, x):
    """

    This function gives you the 2D interpolated result between two points.

    Args:
        a: array with the two x points you want to interpolate between.
        b: array with the two y points you want to interpolate between.
        x: point on the x axis at which we want to interpolate.

    Usage:
        >> y_interpol = interpolation(a, b, x)
    
    """

    output = b[0] + (x - a[0]) * ((b[1] - b[0])/(a[1] - a[0]))
    return output

def find(name, path):
    """

    This function finds the document stated.
    If the document is found it returns the whole path name, if it is not found it returns 'None'.

    Args:
        name: file name (string).
        path: file path (string).

    Usage:
        >> path = find(name, path)
    
    """

    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def MS_turn_off():
    """

    Calculates the age at which low mass stars (0.1-1.3 Msol) leave the Main Sequence

    Args:
        None.

    Usage:
        >> lifetimes = MS_turn_off()
    
    """
    t_sol=10 #Main sequence lifetime of the Sun (in Gyrs)
    m_sol=1 #Mass of Sun in solar masses

    #To find lifetime of a star in the Main Sequence: t \proportional 1/M^(2.5)
    M=MIST_masses()

    lifetime=[]

    for i in M:
        t=t_sol*(1/i**(2.5))
        lifetime.append(t)
        #print('A star with mass ', i, ' M_sol will have a Main Sequence lifetime of ', t, ' Gyrs.')
        
    return lifetime

def Stop_MIST():
    """

    Returns the age at which the MIST tables should stop calculating the evolution of the star.

    Args:
        None.

    Usage:
        >> Stop_MIST = Stop_MIST()
    
    """
    stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,
              20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,
              12317,11000,10000,8500,8000,7500,7000,6500,6000,5500,4500,4000,3500,3000]
        
    return stop_MIST

def crit_age():
    """

    Returns the age at which low mass stars start to behave weirdly.

    Args:
        None.

    Usage:
        >> age_crit = crit_age()
    
    """
    age_crit=[230*1e6, 250*1e6, 270*1e6, 290*1e6]
        
    return age_crit
