'For each mass display the evolution of 1 day rotation period and when their Main Sequence lifetime ends'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/' )

from Final_MIST_tables import common_lib
from Final_MIST_tables import spindown_model
from Final_MIST_tables import load_mist_models

def constrained_Evolution_plot(x, lim_period, Prot_evol):
    """

    Plots the rotation period evolution of all MIST masses (with initial rotation period of 1 day) and the Main Sequence turn-off of all these stars.

    Args:
        x: Main Sequence lifetime of star.
        lim_period: period limits (y-axis) in plot.
        Prot_evol: array that contains the rotation period evolution of the star for each age step.

    Usage:
        >> constrained_Evolution_plot(x, lim_period, Prot_evol)
    
    """
    #Create evolution plot
    fig = plt.figure(figsize=(10, 5))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) 
    plt.title(f'P_rot evolution of {MASSES[i]} $M_\odot$ star with P_rot,i = 1 day', fontsize= 12)
    plt.xlim([10,(lifetime[i]*1000)+10000])
    plt.ylim([ 0.01, 500])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.plot(x, lim_period, color='red') #Plot the age at which the star leaves the main sequence
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2])
    #plt.savefig(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Main_Sequence_Lifetime_{MASSES[i]}.png', bbox_inches='tight')


##### BODY OF CODE STARTS HERE #####

Prot=np.ones(1)*1 #Rotation period 1 day
MASSES=common_lib.MIST_masses() #exact masses
lifetime=common_lib.Stop_MIST() #when to stop the star evolution

#Create rotation period evolution plots
for i in range(len(MASSES)):
    used_M_resampled=MASSES[i] #Take mass from mass list
    M0=np.ones(1)*used_M_resampled #Create an array with that mass value

    Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=Prot, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
    
    lim_period=np.arange(0.01, 500,0.1) #Period limits in plot
    x=np.ones((lim_period.shape))*lifetime[i]#*1000 #Main Sequence lifetime of star

    constrained_Evolution_plot(x, lim_period, Prot_evol)
    
    