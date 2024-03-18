import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import astropy.units as u
import astropy.constants as const
import sys

'For a mass of 0.1-0.3 solar mass derive the rotation period evolution for a synthetic data of initial rotation periods'

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

import common_lib
import load_mist_models
import spindown_model

def plot_Small_M_Prot_evolution(Prot_evol):
    
    """
    Plots the rotation period evolution of a star with a mass = 0.1-0.3M_sol for different initial rotation periods.

    Args:
        Prot_evol: array containing all rotation period values at different stages in the star's life.

    Usage:
        >> plot_Small_M_Prot_evolution(Prot_evol):
    
    """
    fig = plt.figure(figsize=(8, 5))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
    plt.xlim([10,0.7*10**4])
    plt.ylim([ 0.01, 10000])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title(f'Stellar spindown plot for {i}$M_\odot$')
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2])
        
        
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #plt.savefig('/home/farah/Documents/Project/Data/Small_M_Evolution_Prot.png')
    

##### BODY OF CODE STARTS HERE ######

# define initial rotation periods
P0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2, 3, 4, 5., 6, 7, 8, 9, 10.])

mass = [0.16]

for i in mass:
    # define initial stellar masses
    M0 = np.ones((P0.shape))*i
    
    # t0: time [Myr] at which the model should start 
    # tdisc: disk-locking time [Myr]
    Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=P0, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
    
    plot_Small_M_Prot_evolution(Prot_evol)
