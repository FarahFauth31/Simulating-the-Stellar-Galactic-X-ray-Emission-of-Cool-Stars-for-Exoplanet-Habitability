import numpy as np
import astropy.units as u
import astropy.constants as const
import sys
import matplotlib.pyplot as plt

'MIST table plots'

sys.path.append( '/Project/MIST_tables' )

def load_mist_tables(Mstar=1., filepath='/Project/MIST_tables'):
        """
        Load in the MIST tables.
        
        Args:
            Mstar: Stellar masses in units of solar masses
            filepath: Path where the MIST tables are stored
        Usage:
            >> AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist=load_mist_tables(Mstar=mass)
            
        """
        import v2_read_mist_models
        import astropy.units as u

        print(filepath+f'/{Mstar}M_history.data')

        eep = v2_read_mist_models.EEP(filepath+f'/{Mstar}M_history.data', verbose=False)
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Myr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['total_moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs
        MASS_mist = eep.eeps['star_mass']
        log_RADIUS_mist = eep.eeps['log_R']
        RADIUS_mist = 10**log_RADIUS_mist

        return AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist

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

# mass01_04 = np.arange(10, 40, 2) / 100 #0.02 steps
# mass04_09 = np.arange(40, 90, 5) / 100 #0.05 steps
# mass09_11 = np.arange(90, 110, 2) / 100 #0.02 steps
# mass11_13 = np.arange(110, 135, 5) / 100 #0.05 steps
# mass = np.concatenate((mass01_04, mass04_09, mass09_11, mass11_13))

mass = np.arange(10, 20, 2) / 100

n_steps=len(mass)

for i in range(n_steps):
    AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist=load_mist_tables(Mstar=mass[i])
    
    #Tau plot
    plot_tau(AGE_mist, TAU_mist, mass[i])
    
    #MoI plot
    plot_MoI(AGE_mist, MOI_mist, mass)
    
    #Mass plot
    plot_mass(AGE_mist, MASS_mist, mass)
    
    #Radius plot
    plot_radius(AGE_mist, RADIUS_mist, mass)
