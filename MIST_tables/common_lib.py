import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

'Common functions'

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