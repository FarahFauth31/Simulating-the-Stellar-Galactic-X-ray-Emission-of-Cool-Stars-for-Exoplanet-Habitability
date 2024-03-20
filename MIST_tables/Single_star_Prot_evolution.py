'Rotation evolution of a single star'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

import spindown_model
import common_lib

#Make inital rotation period distribution and resample Prots for different star masses
def initial_Prot_dist():
    
    """
    Creates an initial rotation period distribution from which it resamples rotation periods for different star masses.
    It returns the resampled rotation period.

    Args:
        None.

    Usage:
        >> Prot_resampled = initial_Prot_dist()
    
    """   
    #Create an initial rotation period distribution from hPer data
    Full_m=data.Mass #Select the data set you want, in this case the mass (y data)
    Full_Prot=data.Per #Select the data set you want, in this case the rotation periods

    #Create new lists. The mass data is limited between 0.3 < m < 1.2 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
    m=[]
    Prot=[]

    for i in range(len(Full_m)):
        if 0.3 <= Full_m[i] <= 1.2 and 0.01 <= Full_Prot[i] <= 12:
            m.append(Full_m[i])
            Prot.append(Full_Prot[i])
            
    #Resample data limiting mass and rotation period again        
    n_sample=1 #Sample number
    step=1 #One data point at a time
    c=0 #Count

    #Empty resampled data arrays with sample number size
    Prot_resampled=np.zeros(n_sample) 
    M_resampled=np.zeros(n_sample)

    #Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
    while c < n_sample:
        values = np.vstack([Prot, m]) # create 2D array that contains the properties you want to resample
        kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
        
        #Create a resampled data point
        re_Prot = kde_ProtMass.resample(step)[0,:]
        re_M = kde_ProtMass.resample(step)[1,:]
        
        if 0.3 <= re_M <= 1.2 and 0.01 <= re_Prot <= 12:
            M_resampled[c]=re_M
            Prot_resampled[c]=re_Prot
            c+=1
    return Prot_resampled


def plot_Prot_evolution(age, Calculated_Prot, Prot_evol, solar_mass):
    
    """
    Plots the rotation period evolution of a star with a mass=solar_mass and an initial rotation period given by initial_Prot_dist().
    It also plots the calculated rotation period the star has at a given age.

    Args:
        age: star age in Myrs.
        Calculated_Prot: rotation period calculated at age.
        Prot_evol: array containing all rotation period values at different stages in the star's life.
        solar_mass: mass of star in units of M_sol.

    Usage:
        >> plot_Prot_evolution(age, Calculated_Prot, Prot_evol, solar_mass)
    
    """   
    #Plot Age vs Rotation period
    fig = plt.figure(figsize=(15, 10))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
    plt.xlim([10,1*10**4])
    plt.ylim([ 0.01, Calculated_Prot + 50])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.scatter(age,Calculated_Prot, c='black', marker='$\odot$', s=500)
    plt.title(f'Evolution of the rotation period of {solar_mass}$M_\odot$ star', fontsize= 20)
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2], c='orange')
        
        
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)

    plt.savefig('/home/farah/Documents/Project/Data/Star_Evolution_Prot.png')

###### BODY OF CODE STARTS HERE ######

data = pd.read_csv("/home/farah/Documents/Project/Data/hPer_Data.csv") #Read hPer csv file
star_data = pd.read_csv("/media/farah/T7 Shield/GaiaUniverseModel_0000.csv") #Read GUMS csv file

MASSES = common_lib.MIST_masses()

# Resampled Prot data
Prot_resampled = initial_Prot_dist()

#Mass of the object we are studying using star data
n_star=5
solar_mass=1#star_data.mass[n_star-1]
used_mass=common_lib.find_nearest(MASSES, solar_mass)
M0 = np.ones((Prot_resampled.shape))*used_mass

# t0: time [Myr] at which the model should start 
# tdisc: disk-locking time [Myr]
Prot_evol, age_zero, Prot_interp, spl_Prot = spindown_model.spin_down_evol(Prot_init=Prot_resampled, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
#Select the rotation period for the age of the star
age=4600#star_data.age[n_star-1]*(10**3) #Age of star (it is given in Gyr and we pass it to Myr)
Calculated_Prot=spl_Prot(age)

#Plot Age vs Rotation period
plot_Prot_evolution(age, Calculated_Prot, Prot_evol, solar_mass)

print('This star with a mass of ',solar_mass,' solar masses has an age of ',age,' Myrs. The calculated rotation period for that age is ',Calculated_Prot,' days.')

