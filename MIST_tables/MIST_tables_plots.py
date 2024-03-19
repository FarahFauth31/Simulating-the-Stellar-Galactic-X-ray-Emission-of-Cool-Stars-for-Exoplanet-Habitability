'MIST table plots'

import sys
import matplotlib.pyplot as plt

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

import load_mist_models
import common_lib

def plot_tau(x, y, mass):
    """
    Plots the convective turnover time against stellar age.

    Args:
        x: stellar age in units of Myrs
        y: convective turnover time in units of days
        mass: the stellar mass in units of solar masses

    Usage:
        >> plot_tau(AGE_mist, TAU_mist, mass)
    
    """    
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    # plt.ylim([ 0, 500])
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('$\u03C4$ [days]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of $\u03C4$ of {mass}$M_\odot$ star', fontsize= 20)
    plt.scatter(x, y, c='blue')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
def plot_MoI(x, y, mass):
    """
    Plots the moment of inertia against stellar age.

    Args:
        x: stellar age in units of Myrs
        y: moment of inertia in cgs units
        mass: the stellar mass in units of solar masses

    Usage:
        >> plot_MoI(AGE_mist, MOI_mist, mass)
    
    """    
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    plt.ylim([ 0, 2])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Moment of Inertia [$x 10^{54} g*s^2$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of MoI of {mass}$M_\odot$ star', fontsize= 20)
    plt.scatter(x, y/(10**(54)), c='red')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)

def plot_mass(x, y, mass):
    """
    Plots the stellar mass against stellar age.

    Args:
        x: stellar age in units of Myrs
        y: mass in units of solar masses
        mass: the stellar mass in units of solar masses

    Usage:
        >> plot_mass(AGE_mist, MASS_mist, mass)
    
    """    
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    #plt.ylim([ 0, 10000000000])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Mass [$M_\odot$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of mass of {mass}$M_\odot$ star', fontsize= 20)
    plt.scatter(x, y, c='green')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)

def plot_radius(x, y, mass):
    """
    Plots the stellar radius against stellar age.

    Args:
        x: stellar age in units of Myrs
        y: stellar radius in units of solar radii
        mass: the stellar mass in units of solar masses

    Usage:
        >> plot_radius(AGE_mist, RADIUS_mist, mass)
    
    """    
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    plt.ylim([ 0, 2])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Radius [$R_\odot$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of radius of {mass}$M_\odot$ star', fontsize= 20)
    plt.scatter(x, y, c='yellow')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)


####### BODY OF CODE STARTS HERE #######

mass = common_lib.MIST_masses()
#mass = np.arange(98, 100, 2) / 100

n_steps=len(mass)

for i in range(n_steps):
    AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist=load_mist_models.load_mist_tables(Mstar=mass[i])
    
    #Tau plot
    plot_tau(AGE_mist, TAU_mist, mass[i])
    
    #MoI plot
    plot_MoI(AGE_mist, MOI_mist, mass[i])
    
    #Mass plot
    plot_mass(AGE_mist, MASS_mist, mass[i])
    
    #Radius plot
    plot_radius(AGE_mist, RADIUS_mist, mass[i])

