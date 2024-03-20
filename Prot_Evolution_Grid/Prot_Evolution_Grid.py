'Rotation period and age evolution of stars - Grid attempt'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import astropy.units as u
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/' )

from Final_MIST_tables import common_lib
from Final_MIST_tables import spindown_model
from Final_MIST_tables import load_mist_models

def empty_Final_set(MASSES, M_new, MIST_data, p_sample):
    """

    Creates empty array with the desired shape [number of initial rotation periods][2][length of MIST data until MS turn-off].

    Args:
        MASSES: array containing MIST masses.
        M_new: array with a single data point, the mass of the star.
        MIST_data: all the data extracted from the MIST file.
        p_sample: length of the rotation period array (how many rotation periods we are looking at).

    Usage:
        >> Final_set = empty_Final_set(MASSES, M_new, MIST_data, p_sample)
    
    """
    #Age at which the MIST table should stop depending on mass of star
    index=int(np.where(MASSES==M_new[0])[0])
    stop_MIST=common_lib.Stop_MIST()
    max_age=stop_MIST[index]*1e6*u.yr  
    # find index at which we want to stop the rotational evolution
    prot_ind = (MIST_data[0][MIST_data[0] < max_age]).shape[0]
    Final_set=np.zeros((p_sample,2,len(MIST_data[0][:prot_ind]))) #Empty array for final data with rotation period evolution over time
    return Final_set

def scatter_plot_Prot_evolution(MASSES, Final_set):
    """

    Plots a scatter plot to study time steps and rotation period evolution for each initial rotation period.

    Args:
        MASSES: array containing MIST masses.
        Final_set: array that contains the rotation period evolution and all the time steps for each initial rotation period given.

    Usage:
        >> scatter_plot_Prot_evolution(Final_set)
    
    """

    #Plot Age vs Rotation period
    fig = plt.figure(figsize=(15, 10))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
            
    plt.xlim([10,2.5*10**4])
    plt.ylim([ 0.01, 10000])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title(f'Evolution of the rotation period of stars - {MASSES[c]}$M_\odot$', fontsize= 25)
    for el in range(p_sample):
        plt.scatter(Final_set[el][0], Final_set[el][1])
            
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #plt.savefig(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Scatter_Time_evolution_study_{MASSES[c]}.png')


##### BODY OF CODE STARTS HERE #####

MASSES = [1.06]#common_lib.MIST_masses()[26:]
Prot=np.concatenate((np.arange(0.1,1,0.01),np.arange(0.1,12.1,0.1))) #0.01 step from 0.1 to 1 days (since it is logarithmic) and 0.1 steps from 1 to 12 days
              
m_sample=len(MASSES) #How many masses we want to look at
p_sample=len(Prot) #How many rotation periods we want to look at

#Empty data one unit arrays so they can be used by the model
Prot_new=np.zeros(1) 
M_new=np.zeros(1)

#Creating a loop to create the evolutionary data we need by taking each initial period and running the model for each mass
for c in range(m_sample): #For each mass
    M_new[0]=MASSES[c] #Put mass number in the single unit mass array
    MIST_data=load_mist_models.load_mist_tables(Mstar=MASSES[c]) #Look at how many time steps the MIST table has
    Final_set=empty_Final_set(MASSES, M_new, MIST_data, p_sample) #Empty array for final data with rotation period evolution over time

    for i in range(p_sample): #For each initial rotation period
        Prot_new[0]=Prot[i] #Put rotation period number in the single unit Prot array
        # t0: time [Myr] at which the model should start 
        # tdisc: disk-locking time [Myr]
        Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=Prot_new, 
                                                                Mstar_init=M_new, 
                                                                t0=1., tdisc=13.)

        Final_set[i][0]=Prot_evol[0,1] #Evolution of time
        Final_set[i][1]=Prot_evol[0,2] #Evolution of rotation period
    
    #scatter_plot_Prot_evolution(MASSES, Final_set)
    
    with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Grid_docs/{MASSES[c]}M_evolution','wb') as f: pickle.dump(Final_set, f) #Save data for one mass as a pickle file
    f.close()

#print(Final_set) #Print final array with all data for one mass

