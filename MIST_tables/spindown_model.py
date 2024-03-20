'Spindown Model (Variation of Garraffo et al. 2018)'

import numpy as np
import sys
import pickle
import astropy.constants as const
import astropy.units as u

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' ) #You can download MIST tables from https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data

import common_lib
import load_mist_models

MASSES = common_lib.MIST_masses()

#Function that defines cgs constants
def define_cgs_constants():
    
    M_sun = const.M_sun.cgs     # solar mass [gr]
    R_sun = const.R_sun.cgs     # solar radius [cm]
    G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

    return M_sun, R_sun, G

#Function that evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
def spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=4e41, n_min=1., complexity=True):

    """
    Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
    Prot_init: Starting rotation periods
    Mstar_init: Corresponding stellar masses
    t0: Time in Myr at which the rotational evolution should be started
    a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
    """
    from scipy.interpolate import InterpolatedUnivariateSpline
    
    M_sun, R_Sun, G = define_cgs_constants()

    # DEFINE PARAMETERS
    Prot0 = Prot_init*u.d               # initial rotation periods
    Mstar = Mstar_init*u.Msun           # initial stellar masses
    N = len(Prot_init)                  # sample size (number of stars)

    # age_max = 1e10*u.yr                 # oldest cluster's age [million years]
    # tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    # ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC
    
    stop_MIST=common_lib.Stop_MIST() # age at which we want to stop the rotation period evolution for each individual mass [Myrs]

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist = load_mist_models.load_data_spindownmodel(Mstar=Mstar[j].value)
        
        #Age from MIST table at which the code should stop depending on mass of star
        index=int(np.where(MASSES==Mstar_init[j])[0]) #Find index of mass in the MASSES array
        max_age=stop_MIST[index]*1e6*u.yr #Use that index to find the age at which the code should stop depending on mass of star
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))

        # Rotation Periods [d]
        Prot = np.tile(Prot0[j], len(AGE_mist))

        # Angular Momentum [g cm^2/d]
        J = MOI_mist*2.*np.pi/Prot  

        # Initial Angular Velocity [1/d]
        Omega = 2.*np.pi/Prot[0]

        ## modulate the angular momentum loss and evolve Prot
        ## go through each time step in the MIST table starting from t_DL
        for i in range(ind0+1, (AGE_mist[AGE_mist < max_age]).shape[0]-1): 
    
            # define timestep [d]
            dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 

            # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
            Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 

            # Rossby number
            Ro = Prot[i]/TAU_mist[i]

            # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
            n = a/Ro + b*Ro + n_min 
            # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
            
            # In the following, consider only stars with magnetic complexity 1 < n < 8
            if (n < 1) or (complexity == False):
                # Turn this on to ignore complexity: Skumanich evolution
                n = 1
            elif n > 8:
                n = 8 

            # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
            B = 100. #Gauss
            QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless

            # Modulate angular momentum loss
            Jdot =  QJ*Jdot_dip
            Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
            Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
            Omega = Omega - dt*Jdot/MOI_mist[i]
            J[i] = MOI_mist[i]* Omega

        # stellar mass [Msun], time [Myr], evolved Prot [d]
        Prot[Prot < 0.] = np.nan
        Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])
        
        ### For masses m=0.1-0.16 [Solar masses], the rotation period evolution behaves weirdly with initial rotation periods P_rot,i < 1.4 days. That's why we smooth the plots
        if Mstar[j].value < 0.18 and Prot_init[j] < 1.4: # for mass < 0.18 M_sol and Prot < 1.4 days, we smooth the plot 
            c1=AGE_mist.value
            crit_age=common_lib.crit_age()[int(np.where(MASSES==Mstar[j].value)[0])] #Age at which stars start to behave weirdly
            d1=common_lib.find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 1.4 days (Create that file running the code for the Prot_init=np.array([1.4]) and Mstar_init=np.array([0.1-0.16]) and saving the Prot ouput as a pickle file)
            with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Prot_i_1.4d_evolution_{Mstar[j].value}M','rb') as f: file = pickle.load(f)
            Prot_i_14=file*u.d            
            Prot[ind:]=Prot_i_14[ind:] #Smooth plot by using the evolution of P_rot,i = 1.4 days after the critical age
            f.close()

        #interpolate all evolved Prots along the same age-array for the direct comparison
        #(as different masses were evolved for different long times in the MIST tables)
        if j==0:
            # define a common age-array for the interpolation
            age_zero = AGE_mist/1e6
        spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
        Prot_interp.append(spl_Prot(age_zero))
        
    Prot_evol = np.array(Prot_evol)   
    Prot_interp = np.array(Prot_interp)
    
    return Prot_evol, age_zero, Prot_interp, spl_Prot


#Function that evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
def break_v_spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=4e41, n_min=1., complexity=True, filter_breakup=True):

    """
    Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
    Prot_init: Starting rotation periods
    Mstar_init: Corresponding stellar masses
    t0: Time in Myr at which the rotational evolution should be started
    a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
    """
    import astropy.units as u
    from scipy.interpolate import InterpolatedUnivariateSpline
    #from scipy.signal import savgol_filter
    
    M_sun, R_Sun, G = define_cgs_constants()

    # DEFINE PARAMETERS
    Prot0 = Prot_init*u.d               # initial rotation periods
    Mstar = Mstar_init*u.Msun           # initial stellar masses
    N = len(Prot_init)                  # sample size (number of stars)

    # age_max = 1e10*u.yr                 # oldest cluster's age [million years]
    # tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    # ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC
    
    stop_MIST=common_lib.Stop_MIST() # age at which we want to stop the rotation period evolution for each individual mass [Myrs]

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist = load_mist_models.load_mist_tables(Mstar=Mstar[j].value)
        
        #Age from MIST table at which the code should stop depending on mass of star
        index=int(np.where(MASSES==Mstar_init[j])[0]) #Find index of mass in the MASSES array
        max_age=stop_MIST[index]*1e6*u.yr #Use that index to find the age at which the code should stop depending on mass of star
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))

        # Rotation Periods [d]
        Prot = np.tile(Prot0[j], len(AGE_mist))

        # Angular Momentum [g cm^2/d]
        J = MOI_mist*2.*np.pi/Prot  

        # Initial Angular Velocity [1/d]
        Omega = 2.*np.pi/Prot[0]
        
        # Breakup rotation period [d]
        Omega_breakup = np.sqrt(G*Mstar[j].cgs/((RADIUS_mist*u.Rsun).cgs**3.))
        Prot_breakup = (2.*np.pi/Omega_breakup).to(u.d)
        Prot_fastrotators = np.zeros(len(Prot))*u.d

        ## modulate the angular momentum loss and evolve Prot
        ## go through each time step in the MIST table starting from t_DL
        for i in range(ind0+1, (AGE_mist[AGE_mist < max_age]).shape[0]-1): 
    
            # define timestep [d]
            dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 

            # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
            Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 

            # Rossby number
            Ro = Prot[i]/TAU_mist[i]

            # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
            n = a/Ro + b*Ro + n_min 
            # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
            
            # In the following, consider only stars with magnetic complexity 1 < n < 8
            if (n < 1) or (complexity == False):
                # Turn this on to ignore complexity: Skumanich evolution
                n = 1
            elif n > 8:
                n = 8 

            # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
            B = 100. #Gauss
            QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless

            # Modulate angular momentum loss
            Jdot =  QJ*Jdot_dip
            Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
            Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
            Omega = Omega - dt*Jdot/MOI_mist[i]
            J[i] = MOI_mist[i]* Omega
        
        time_filter = (AGE_mist > t0) & (AGE_mist < max_age)
        if (filter_breakup == True) and (np.any(Prot[time_filter] <= Prot_breakup[time_filter])):
            ind_Prot_breakup = np.array(np.where(Prot[time_filter] <= Prot_breakup[time_filter]))+len(Prot)-len(Prot[time_filter])
            Prot_fastrotators = Prot
            Prot = np.zeros(len(Prot))*u.d
            Prot[:] = np.nan

        # stellar mass [Msun], time [Myr], evolved Prot [d]
        Prot[Prot < 0.] = np.nan
        Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot, Prot_breakup])
        
        ### For masses m=0.1-0.16 [Solar masses], the rotation period evolution behaves weirdly with initial rotation periods P_rot,i < 1.4 days. That's why we smooth the plots
        if Mstar[j].value < 0.18 and Prot_init[j] < 1.4: # for mass < 0.18 M_sol and Prot < 1.4 days, we smooth the plot 
            c1=AGE_mist.value
            crit_age=common_lib.crit_age()[int(np.where(MASSES==Mstar[j].value)[0])] #Age at which stars start to behave weirdly
            d1=common_lib.find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 1.4 days (Create that file running the code for the Prot_init=np.array([1.4]) and Mstar_init=np.array([0.1-0.16]) and saving the Prot ouput as a pickle file)
            with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Prot_i_1.4d_evolution_{Mstar[j].value}M','rb') as f: file = pickle.load(f)
            Prot_i_14=file*u.d            
            Prot[ind:]=Prot_i_14[ind:] #Smooth plot by using the evolution of P_rot,i = 1.4 days after the critical age
            f.close()

        #interpolate all evolved Prots along the same age-array for the direct comparison
        #(as different masses were evolved for different long times in the MIST tables)
        if j==0:
            # define a common age-array for the interpolation
            age_zero = AGE_mist/1e6
        spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
        Prot_interp.append(spl_Prot(age_zero))
        
    Prot_evol = np.array(Prot_evol)   
    Prot_interp = np.array(Prot_interp)
    
    return Prot_evol, age_zero, Prot_interp, spl_Prot
