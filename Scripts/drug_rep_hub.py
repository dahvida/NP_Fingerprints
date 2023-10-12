from sim_search import *
from utils import *
import pandas as pd
import numpy as np
import os
from multiprocessing import Process, Pool
import argparse

######################################################################

def main():
    
    #get names of all FPs available in ../FPs
    names = os.listdir("../FPs/drug_rep_hub")
    fp_names = [x[:-4] for x in names]
    fp_paths = ["../FPs/drug_rep_hub/" + x for x in names]
    
    
    
    #compute number of pairwise similarity calculations and
    #preallocate array of correct size
    total_pairs = int((6776**2 - 6776) / 2 )
    pair_matrix = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    stats_matrix = np.zeros((22, len(fp_paths)),
                             dtype = np.float32)
    
    #create paths to save files
    path = "../Results/drug_rep_hub/"
        
    #loop over each FP in ../FPs
    for j in range(len(fp_paths)):
        
        #load FP
        print(f"[drug_rep_hub]: Currently processing {fp_paths[j]}...")
        fp = pickle_load(fp_paths[j], verbose=False)

        #depending on whether it is minhashed or not, run search
        if fp_names[j] != "mhfp" and fp_names[j] != "map4":
            sims = sim_search(fp)
        else:
            sims = sim_search(fp, sim_metric="custom")
        
        #get statistics of similarity distribution
        stats = eval_sim(sims)
        
        #save i-th statistics and i-th pairwise sims in arrays
        stats_matrix[:,j] = stats
        pair_matrix[j,:] = sims
    
    #calculate correlations between FPs
    print(f"[drug_rep_hub]: Calculating correlation matrix...")
    correlation_matrix = np.corrcoef(pair_matrix)
    
    #save everything
    save_corr_df(correlation_matrix, fp_names,
                     path = path + "corr.csv",
                     verbose = False)
    save_sim_df(stats_matrix, fp_names,
                    path = path + "sim.csv",
                    verbose = False)


if __name__ == "__main__":
    main()

