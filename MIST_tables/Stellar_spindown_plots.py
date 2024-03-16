import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import astropy.units as u
import astropy.constants as const
import sys

'For a random mass sample create synthetic rotation period data and plot the rotation evolution for each mass'

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

import common_lib
import load_mist_models
import spindown_model

#Create initial rotation period data with hPer data
def random_star_sample(hPer_data):
    
    """
    Creates a kde gaussian distribution with cluster data and then it resamples mass and rotation period values for a synthetic population of stars.

    Args:
        hPer_data: cluster star data.

    Usage:
        >> Prot_resampled, M_resampled = random_star_sample(hPer_data)
    
    """
    m=hPer_data.Mass #Select the data set you want, in this case the mass (y data)
    Prot=hPer_data.Per #Select the data set you want, in this case the rotation periods
    l_Prot=np.log10(Prot) #Make the logarithm of selected data (x data)

    n_sample=4
    values_Mass = np.vstack([Prot, m]) # create 2D array that contains the properties you want to resample
    Mass_new = np.linspace(0.3,1.2,n_sample) # define a new array in which you want to resample
    kde_ProtMass = stats.gaussian_kde(values_Mass) # calculate 2D KDE of Prot-vs-Mstar

    Prot_resampled = kde_ProtMass.resample(n_sample)[0,:]
    M_resampled = kde_ProtMass.resample(n_sample)[1,:]
    
    return Prot_resampled, M_resampled

def plot_Prot_evolution(Prot_evol, M_resampled, i):
    
    """
    Plots the rotation period evolution of a star with a mass=solar_mass for different initial rotation periods given by random_star_sample(hPer_data).

    Args:
        Prot_evol: array containing all rotation period values at different stages in the star's life.
        M_resampled: resampled mass of star in units of M_sol.
        i: index of star in sample.

    Usage:
        >> plot_Prot_evolution(Prot_evol, M_resampled, i)
    
    """
    fig = plt.figure(figsize=(8, 5))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
    
    plt.xlim([10,0.7*10**4])
    plt.ylim([ 0.01, 200])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title(f'Stellar spindown plot ({M_resampled[i]}$M_\odot$)')
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2])
    
    
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    #plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)

###### BODY OF CODE STARTS HERE #####

#Create mass array to find nearest neighbour
MASSES = common_lib.MIST_masses() #all the MIST masses we have

hPer_data = pd.read_csv("/home/farah/Documents/Project/Data/hPer_Data.csv") #Read csv file

Prot_resampled, M_resampled = random_star_sample(hPer_data)

for i in range(len(M_resampled)):
    used_M_resampled=common_lib.find_nearest(MASSES, M_resampled[i])
    M0=np.ones((Prot_resampled.shape))*used_M_resampled

    # t0: time [Myr] at which the model should start 
    # tdisc: disk-locking time [Myr]
    Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=Prot_resampled, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
    
    plot_Prot_evolution(Prot_evol, M_resampled, i)
