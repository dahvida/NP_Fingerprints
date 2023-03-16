"""
General utility functions used throughout the repository
"""

import pandas as pd
import numpy as np
import pickle
#suppress warnings due to 0 division when dealing with unlucky data samples
#(i.e. if there are no compounds of class marine in the i-th sample)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

###########################################################################

def get_statistics(
        df: pd.DataFrame
        ) -> pd.DataFrame:
    """Collects class and scaffold stats for a database slice

    Args:
        df:     dataframe to process for the analysis

    Returns:
        Dataframe with number of compounds, dataset %, number of
        scaffolds and scaffold diversity % for each organism origin
        class
    """

    #prealloc data with right size
    data = np.zeros((6,4))
    
    #make organism labels
    organisms = [
                "all",
                "plant",
                "bacteria",
                "fungi",
                "animal",
                "marine"]
    
    #make column labels
    columns = ["Number of compounds",
               "Dataset %",
               "Number of scaffolds",
               "Scaffold diversity %"]
    
    #collect stats when considering all classes together
    data[0,0] = len(df)
    data[0,1] = 100.0

    scaffs = list(df["murko_framework"])
    data[0,2] = len(set(scaffs))
    data[0,3] = data[0,2] * 100 / data[0,0]
    
    #catch exception if class is absent
    try:
        #loop over all organisms
        for i in range(1, len(organisms)):
            #slice according to class
            temp_df = df.loc[df["organism"] == organisms[i]]
            
            #collect stats
            data[i,0] = len(temp_df)
            data[i,1] = data[i,0] * 100 / data[0,0]
            scaffs = list(temp_df["murko_framework"])
            data[i,2] = len(set(scaffs))
            data[i,3] = data[i,2] * 100 / data[i,0]
    except:
        pass

    #assemble df
    output = pd.DataFrame(
            data = data,
            index = organisms,
            columns = columns
            )
    
    return output

#--------------------------------------------------------------------------#

def pickle_save(array,
                path,
                verbose = True):
    """
    Saves numpy array to path as .pkl file
    """
    with open(path, "wb") as output_file:
        pickle.dump(array, output_file)
    if verbose is True:
        print(f"[utils]: File saved at {path}")

#--------------------------------------------------------------------------#

def pickle_load(path,
                verbose = True):
    """
    Loads .pkl file from path
    """
    with open(path, "rb") as file:
        output = pickle.load(file)
    if verbose is True:
        print(f"[utils]: Loaded {path}")
    return output








