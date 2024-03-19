'For a mass of 1 solar mass derive the rotation period evolution for a synthetic data of initial rotation periods'

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

import common_lib
import spindown_model

def plot_1M_Prot_evolution(Prot_evol):
    
    """
    Plots the rotation period evolution of a star with a mass = 1M_sol for different initial rotation periods.

    Args:
        Prot_evol: array containing all rotation period values at different stages in the star's life.

    Usage:
        >> plot_1M_Prot_evolution(Prot_evol)
    
    """
    colors = plt.cm.spring(np.linspace(0,1,len(P0)))
    fig = plt.figure(figsize=(8, 5))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
    plt.xlim([10,0.7*10**4])
    plt.ylim([ 0.01, 200])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title('Stellar spindown plot for 1$M_\odot$')
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2],c=colors[el])

    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #plt.savefig('/home/farah/Documents/Project/Data/1M_Evolution_Prot.png')

###### BODY OF CODE STARTS HERE #####

#Create mass array to find nearest neighbour
MASSES = common_lib.MIST_masses() #all the MIST masses we have

# define initial rotation periods
P0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2, 3, 4, 5., 6, 7, 8, 9, 10.])

# define initial stellar masses
M0 = np.ones((P0.shape))*1

# t0: time [Myr] at which the model should start 
# tdisc: disk-locking time [Myr]
Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=P0, 
                                                    Mstar_init=M0, 
                                                    t0=1., tdisc=13.)
plot_1M_Prot_evolution(Prot_evol)
