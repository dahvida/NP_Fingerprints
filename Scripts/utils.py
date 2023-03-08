import pandas as pd
import numpy as np

###########################################################################

def get_statistics(df):

    data = np.zeros((6,4))
    
    organisms = [
                "all",
                "plant",
                "bacteria",
                "fungi",
                "animal",
                "marine"]

    columns = ["Number of compounds",
               "Dataset %",
               "Number of scaffolds",
               "Scaffold diversity %"]
    
    data[0,0] = len(df)
    data[0,1] = 100.0

    scaffs = list(df["murko_framework"])
    data[0,2] = len(set(scaffs))
    data[0,3] = data[0,2] * 100 / data[0,0]
    
    for i in range(1, len(organisms)):
        temp_df = df.loc[df["organism"] == organisms[i]]
        data[i,0] = len(temp_df)
        data[i,1] = data[i,0] * 100 / data[0,0]
        scaffs = list(temp_df["murko_framework"])
        data[i,2] = len(set(scaffs))
        data[i,3] = data[i,2] * 100 / data[i,0]

    output = pd.DataFrame(
            data = data,
            index = organisms,
            columns = columns
            )
    
    return output














