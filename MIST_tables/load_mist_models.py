'Load MIST models'

import astropy.units as u
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/' )

def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/Farah_EEPfiles'):
        """
        Load in the MIST tables.
        
        Args:
            Mstar: Stellar masses in units of solar masses
            filepath: Path where the MIST tables are stored
        Usage:
            >> AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist=load_mist_tables(Mstar=mass)
            
        """
        import v2_read_mist_models
    
        print(filepath+f'/{Mstar}M_history.data.eep')

        eep = v2_read_mist_models.EEP(filepath+f'/{Mstar}M_history.data.eep', verbose=False)
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Myr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_turnover_time_l_hybrid']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['total_moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs
        MASS_mist = eep.eeps['star_mass']
        log_RADIUS_mist = eep.eeps['log_R']
        RADIUS_mist = 10**log_RADIUS_mist

        return AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist
    
def load_data_spindownmodel(Mstar=1., filepath='/home/farah/Documents/Redo_Project_Cfa/Final_MIST_tables/Farah_EEPfiles'):
        """
        Load in the MIST tables.
        Mstar: Stellar masses in units of solar masses
        filepath: Path where the MIST tables are stored
        """
        import v2_read_mist_models

        print(filepath+f'/{Mstar}M_history.data.eep')

        eep = v2_read_mist_models.EEP(filepath+f'/{Mstar}M_history.data.eep', verbose=False)
        AGE_mist = eep.eeps['star_age']*u.yr # stellar age in years
        TAU_mist = (eep.eeps['conv_turnover_time_l_hybrid']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['total_moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs

        return AGE_mist, TAU_mist, MOI_mist