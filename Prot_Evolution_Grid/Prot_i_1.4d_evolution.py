'Create document that contains the rotation period evolution for small mass starts with Prot_i = 3 days'
'(this document will help smooth Prot evolution tracks for very low mass stars)'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys
import astropy.units as u
import pickle

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/' )

from Final_MIST_tables import common_lib
from Final_MIST_tables import spindown_model
from Final_MIST_tables import load_mist_models

def scatter_plot_Prot_evolution(Final_set):
    """

    Plots a scatter plot to study time steps and rotation period evolution for each initial rotation period.

    Args:
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
    
    #plt.savefig(f'/home/farah/Documents/Project/Data/Evolution_grid/Scatter_Time_evolution_study_{name_mass[c]}.png')
    
def check_pickle_file(MASSES, Final_set):
    """

    Opens pickle file and checks that the final array is equal to the information stored in the file.

    Args:
        MASSES: mass array.
        Final_set: array that contains the rotation period evolution and all the time steps for each initial rotation period given.

    Usage:
        >> check_pickle_file(MASSES, Final_set)
    
    """

    with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Prot_i_1.4d_evolution_{MASSES[c]}M','rb') as f: yes = pickle.load(f)
        
    b=np.array_equal(Final_set,yes) #Sanity check
    print(b)

    
##### BODY OF CODE STARTS HERE #####

MASSES = [0.1,0.12,0.14,0.16] #array with small masses
RotationP=[1.4] #minimum Prot at which there are no abrupt jumps in Prot
              
m_sample=len(MASSES) #How many masses we want to look at
p_sample=len(RotationP) #How many rotation periods we want to look at

#Empty data one unit arrays so they can be used by the model
Prot_new=np.zeros(1) 
M_new=np.zeros(1)

#Creating a loop to create the evolutionary data we need by taking each initial period and running the model for each mass
for c in range(m_sample): #For each mass
    M_new[0]=MASSES[c] #Put mass number in the single unit mass array
    a=load_mist_models.load_mist_tables(Mstar=MASSES[c]) #Look at how many time steps the MIST table has
    #Age at which the MIST table should stop depending on mass of star
    index=int(np.where(MASSES==M_new[0])[0])
    stop_MIST=common_lib.Stop_MIST()
    max_age=stop_MIST[index]*1e6*u.yr  
    # find index at which we want to stop the rotational evolution
    prot_ind = (a[0][a[0] < max_age]).shape[0]
    Final_set=np.zeros((p_sample,2,len(a[0][:prot_ind]))) #Empty array with final data with rotation period evolution over time

    for i in range(p_sample): #For each initial rotation period
        Prot_new[0]=RotationP[i] #Put rotation period number in the single unit Prot array
        # t0: time [Myr] at which the model should start 
        # tdisc: disk-locking time [Myr]
        Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=Prot_new, 
                                                                Mstar_init=M_new, 
                                                                t0=1., tdisc=13.)

        Final_set[i][0]=Prot_evol[0,1] #Evolution of time
        Final_set[i][1]=Prot_evol[0,2] #Evolution of rotation period
        
    scatter_plot_Prot_evolution(Final_set)
    
    with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Prot_i_1.4d_evolution_{MASSES[c]}M','wb') as f: pickle.dump(Final_set[0][1], f) #Save data for one mass as a pickle file
    #print(Final_set[0][1])
    #check_pickle_file(MASSES, Final_set[0][1])

