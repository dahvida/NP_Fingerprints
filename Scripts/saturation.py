from sim_search import *
from utils import *
import pandas as pd
import numpy as np
import os
from multiprocessing import Process, Pool
import argparse

###############################################################################

parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--drugrephub', action="store_true",
                    help="Whether to use DrugRepHub for the analysis")

parser.add_argument('--coconut', action="store_true",
                    help="Whether to use COCONUT for the analysis")

args = parser.parse_args()

###############################################################################

def get_means(paths, dataset):
    
    fp_names = [x[:-4] for x in paths]
    fp_paths = ["../FPs/" + dataset + "/" + x for x in paths]
    
    means = np.zeros((len(paths), 1))
    
    for i in range(len(fp_paths)):
        fp = pickle_load(fp_paths[i], verbose=False)
        count = np.count_nonzero(fp, axis=1)
        mean_bit = count / fp.shape[1]
        mean_sample = np.mean(mean_bit, axis=0)
        means[i] = mean_sample
    
    output = pd.DataFrame(means, index=fp_names)
    
    return output

#-----------------------------------------------------------------------------#

def drugrephub():
    
    names = os.listdir("../FPs/drug_rep_hub")
    
    df = get_means(names, "drug_rep_hub")
    
    df.to_csv("../Results/saturation/drugrephub.csv")

#-----------------------------------------------------------------------------#
    
def coconut():
    
    names = os.listdir("../FPs/coconut")
    
    df = get_means(names, "coconut")
    
    df.to_csv("../Results/saturation/coconut.csv")
    
#-----------------------------------------------------------------------------#
    
if __name__ == "__main__":
    if args.drugrephub is True:
        drugrephub()
    if args.coconut is True:
        coconut()