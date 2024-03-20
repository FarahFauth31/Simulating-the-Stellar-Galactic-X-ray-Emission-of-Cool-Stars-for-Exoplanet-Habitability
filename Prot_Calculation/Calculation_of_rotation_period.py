'Calculate Prot of GUMS stars for all documents - More efficient'

import pandas as pd
import numpy as np
from scipy import stats
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline
import sys

sys.path.append( '/home/farah/Documents/Redo_Project_Cfa/' )

from Final_MIST_tables import common_lib
from Final_MIST_tables import spindown_model
from Final_MIST_tables import load_mist_models
        
def kde(file_directory):
    """
    
    This function creates a 2D Gaussian kde distribution from cluster data.
    Returns the Prot-Mass 2D Gaussian distribution.

    Args:
        file_directory: Path of the cluster data document.

    Usage:
        >> kde_ProtMass = kde(file_directory)
    
    """
    #Open csv file containing initial rotation period distribution data
    data = pd.read_csv(file_directory) #Read hPer csv file
            
    df = pd.DataFrame(data)
            
    #Create new lists. The mass data is limited between 0.1 < m < 1.25 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
    red_data = df[(df['Mass']>=0.1) & (df['Mass']<=1.25) & (df['Per']>=0.1) & (df['Per']<=12)]
            
    #Create an initial rotation period distribution from hPer data
    Full_m=red_data.Mass #Select the data set you want, in this case the mass (y data)
    Full_Prot=red_data.Per #Select the data set you want, in this case the rotation periods
    values = np.vstack([Full_Prot, Full_m]) # create 2D array that contains the properties you want to resample
    kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
    
    return kde_ProtMass
    
        
def initial_Prot(kde_ProtMass):
    """
    
    This function resamples a initial rotation period value from a 2D Gaussian kde distribution.
    Returns the resampled rotation period value.

    Args:
       kde_ProtMass: Prot-Mass 2D Gaussian distribution.

    Usage:
        >> re_Prot = initial_Prot(kde_ProtMass)
    
    """
    #Resample data limiting mass and rotation period again        
    step=1 #One data point at a time
    c=0 #Count

    #Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
    while c < 1:
        #Create a resampled data point
        re_Prot = kde_ProtMass.resample(step)[0,:]
        re_M = kde_ProtMass.resample(step)[1,:]
        
        if 0.1 <= re_M <= 1.25 and 0.1 <= re_Prot <= 12:
            c+=1
    
    return re_Prot

def open_mass_files():
    """
    
    This function opens all mass files from the Rotation Period Evolution Grid and stores their information in a main array for easy access.
    Returns the main array with all loaded mass documents.

    Args:
       None.

    Usage:
        >> main_array = open_mass_files()
    
    """ 
    #Open grid of the masses
    MASSES = common_lib.MIST_masses()
    main_array=[]
    for i in MASSES:
        with open(f'/home/farah/Documents/Redo_Project_Cfa/Prot_Evolution_Grid/Grid_docs/{i}M_evolution','rb') as f: doc = pickle.load(f)
        main_array.append(doc)
        f.close()
    
    return main_array
    
    
def calculate_prot(Proti_file_directory, GUMS_file_directory):
    """
    
    This function calculates the rotation period of stars in a document from basic stellar parameters and saves original csv file with an extra column for calculated Prots.

    Args:
       Proti_file_directory: File directory of document containing young cluster data aimed for the resampling of initial rotation periods.
       GUMS_file_directory: Path of file with GUMS data or stellar data.

    Usage:
        >> calculate_prot(Proti_file_directory, GUMS_file_directory)
    
    """ 

    'INPUT'
    
    main_array=open_mass_files()
    kde_ProtMass = kde(Proti_file_directory)
    #Create array with initial rotation periods present in the grid
    Initial_Prots=np.arange(0.1,12.1,0.1) #0.1 steps
    #Create mass array with all the MIST masses we have
    MASSES = [0.1,0.15, 0.2, 0.25, 0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15, 1.2, 1.25] #all the MIST masses we have
                    
    #Open csv file containing star data
    mrdata = pd.read_csv(GUMS_file_directory) #Read csv file
    
    'BODY OF THE CODE'
    
    #Create empty lists that we will use to create a new file with the calculated data
    PROT=[]
    
    #Do this loop for each star in file
    n_star=10#len(mrdata['ra']) #Number of stars in file we want to evaluate
                    
    for i in range(n_star):
                                   
        Prot_i=initial_Prot(kde_ProtMass)[0] #Initial rotation period
        Grid_Prot_i1,Grid_Prot_i2=common_lib.find_2_nearest(Initial_Prots, Prot_i) #Nearest initial rotation periods in grid
        index1=int(np.where(Initial_Prots==Grid_Prot_i1)[0]) #Index of 1st nearest initial rotation period
        index2=int(np.where(Initial_Prots==Grid_Prot_i2)[0]) #Index of 2nd nearest initial rotation period
        mass=mrdata.mass[i] #Mass in solar masses
        age=mrdata.age[i]*1000 #uniform_ages[i]*1000 #Age in Myrs
                                
        Mass_1,Mass_2=common_lib.find_2_nearest(MASSES, mass) #Find 2 nearest masses from MIST tables
                        
        main_index1=int(np.where(MASSES==Mass_1)[0]) #Index of 1st nearest initial rotation period
        main_index2=int(np.where(MASSES==Mass_2)[0]) #Index of 2nd nearest initial rotation period
                        
        arrayname1=main_array[main_index1]       
        arrayname2=main_array[main_index2]                        
                        
        #Interpolate rotation period throughout whole life of star for the two MIST tables
        spl_Prot_1 = InterpolatedUnivariateSpline(arrayname1[index1][0], arrayname1[index1][1], ext=0) #Interpolated line rotation period of Mass_1 and P_rot,i_1
        spl_Prot_2 = InterpolatedUnivariateSpline(arrayname2[index1][0], arrayname2[index1][1], ext=0) #Interpolated line rotation period of Mass_2 and P_rot,i_1
        spl_Prot_3 = InterpolatedUnivariateSpline(arrayname1[index2][0], arrayname1[index2][1], ext=0) #Interpolated line rotation period of Mass_1 and P_rot,i_2
        spl_Prot_4 = InterpolatedUnivariateSpline(arrayname2[index2][0], arrayname2[index2][1], ext=0) #Interpolated line rotation period of Mass_2 and P_rot,i_2
                            
        #Calculate rotation period at specific age of star for the two MIST tables for different mass and different initial rotation period
        interp_Prot_1=float(spl_Prot_1(age)) #Interpolated rotation period of Mass_1 and P_rot,i_1
        interp_Prot_2=float(spl_Prot_2(age)) #Interpolated rotation period of Mass_2 and P_rot,i_1
        interp_Prot_3=float(spl_Prot_3(age)) #Interpolated rotation period of Mass_1 and P_rot,i_2
        interp_Prot_4=float(spl_Prot_4(age)) #Interpolated rotation period of Mass_2 and P_rot,i_2
                        
        if np.isnan(interp_Prot_1) or np.isnan(interp_Prot_2) or np.isnan(interp_Prot_3) or np.isnan(interp_Prot_4):
            Final_Prot = np.nan
        else:
            #Create arrays with the two initial rotation periods we are looking at and the rotation periods calculated for each mass
            two_Prot_i=[Grid_Prot_i1,Grid_Prot_i2] #Initial rotation periods
            diff_Prot1=[interp_Prot_1,interp_Prot_3] #For Mass_1
            diff_Prot2=[interp_Prot_2,interp_Prot_4] #For Mass_2
                            
            #Calculate rotation period for each mass for Prot_i by interpolating between the two calculated rotation periods for the two nearest initial rotation periods in grid
            Med_Prot1 = common_lib.interpolation(two_Prot_i, diff_Prot1, Prot_i)
            Med_Prot2 = common_lib.interpolation(two_Prot_i, diff_Prot2, Prot_i)
                            
            #Create lists with the two nearest masses and the calculated rotation periods
            two_masses=[Mass_1,Mass_2]
            two_periods=[Med_Prot1,Med_Prot2]
                            
            #Calculate the final rotation period by interpolating between the two previous results
            Final_Prot = common_lib.interpolation(two_masses, two_periods, mass)
                                
        PROT.append(Final_Prot)
    
    #dictionary = {'Prot': PROT}  
    #dataframe = pd.DataFrame(dictionary) 
    #mrdata['Prot'] = dataframe
    #mrdata.to_csv(GUMS_file_directory, index=False)
    #print(PROT)


##### BODY OF CODE STARTS HERE #####

Proti_file_directory = "/home/farah/Documents/Project/Data/hPer_Data.csv"
GUMS_file_directory = "/media/farah/T7 Shield/GaiaUniverseModel_0000.csv"
calculate_prot(Proti_file_directory, GUMS_file_directory)
