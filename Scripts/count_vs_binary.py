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
                    
parser.add_argument("--n_cores",
                    default = 10,
                    type = int,
                    help = "number of cores to use for the analysis")
                    
parser.add_argument("--n_replicates",
                    default = 50,
                    type = int,
                    help = "number of replicates to use for the analysis")
                    
parser.add_argument("--sample_size",
                    default = 10000,
                    type = int,
                    help = "number of samples to draw for each batch")
                    
parser.add_argument('--random_seed',
                    default = 42,
                    type = int,
                    help = "random seed for reproducibility")

args = parser.parse_args()

###############################################################################

def drugrephub():
    
    #get names of all FPs available in ../FPs
    names = ["ap.pkl", "tt.pkl", "avalon.pkl"]
    fp_names = ["ap", "tt", "avalon", "ap_c", "tt_c", "avalon_c"]
    fp_paths = ["../FPs/drug_rep_hub/" + x for x in names]
    
    #compute number of pairwise similarity calculations and
    #preallocate array of correct size
    total_pairs = int((6776**2 - 6776) / 2 )
    pair_matrix_c = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    pair_matrix_b = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    
    #create paths to save files
    path = "../Results/count_vs_binary/"
        
    #loop over each FP in ../FPs
    for j in range(len(fp_paths)):
        
        #load FP
        print(f"[drug_rep_hub]: Currently processing {fp_paths[j]}...")
        fp = pickle_load(fp_paths[j], verbose=False)

        #depending on whether it is minhashed or not, run search
        sims_b = sim_search(fp)
        sims_c = sim_search(fp, force_binary=False)
        
        #save i-th statistics and i-th pairwise sims in arrays
        pair_matrix_b[j,:] = sims_b
        pair_matrix_c[j,:] = sims_c

    #concat matrices
    pair_matrix = np.concatenate((pair_matrix_b, pair_matrix_c), axis=0)
    
    #calculate correlations between FPs
    print("[drug_rep_hub]: Calculating correlation matrix...")
    correlation_matrix = np.corrcoef(pair_matrix)
        
    #save everything
    save_square_df(correlation_matrix, fp_names,
                         path = path + "c_v_b_drug.csv",
                         verbose = False)
    
#-----------------------------------------------------------------------------#

def coconut(index, iteration_n):
    
    #get relevant index according to iteration number
    print(f"[coconut]: Running iteration {iteration_n}")
    index = index[:,iteration_n]
    n_compounds = index.shape[0]
    
    #get names of all FPs available in ../FPs
    names = ["ap.pkl", "tt.pkl", "avalon.pkl"]
    fp_names = ["ap", "tt", "avalon", "ap_c", "tt_c", "avalon_c"]
    fp_paths = ["../FPs/coconut/" + x for x in names]
    
    #compute number of pairwise similarity calculations and
    #preallocate array of correct size
    total_pairs = int((n_compounds**2 - n_compounds) / 2 )
    pair_matrix_c = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    pair_matrix_b = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    
    #create paths to save files
    path = "../Results/count_vs_binary/"
        
    #loop over each FP in ../FPs
    for j in range(len(fp_paths)):
        
        #load FP
        fp = pickle_load(fp_paths[j], verbose=False)
        
        #overwrite FP with its slice according to index (important to 
        #limit memory usage to the bare minimum we need for i-th run)
        fp = fp[index,:]

        #depending on whether it is minhashed or not, run search
        sims_b = sim_search(fp)
        sims_c = sim_search(fp, force_binary=False)
        
        #save i-th statistics and i-th pairwise sims in arrays
        pair_matrix_b[j,:] = sims_b
        pair_matrix_c[j,:] = sims_c

    #concat matrices
    pair_matrix = np.concatenate((pair_matrix_b, pair_matrix_c), axis=0)
    
    #calculate correlations between FPs
    print("[coconut]: Calculating correlation matrix...")
    correlation_matrix = np.corrcoef(pair_matrix)
        
    #save everything
    save_square_df(correlation_matrix, fp_names,
                         path = path + str(iteration_n) + ".csv",
                         verbose = False)

#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    if args.drugrephub is True:
        drugrephub()
        
    if args.coconut is True:
    	#read params from args
        seed = args.random_seed
        sample_size = args.sample_size
        n_replicates = args.n_replicates
        n_cores = args.n_cores
    	
        #set seed for index generation
        print(f"[coconut]: Beginning analysis")
        np.random.seed(seed)

        #generate matrix to sample from the cleaned dataset for the
        #correlation analysis
        index = np.random.randint(low = 0,
                                  high = 129869,
                                  size = (sample_size, n_replicates))
        
        #save as pkl in ../Results
        pickle_save(index, "../Results/index.pkl", verbose = False)
        
        #create pool with n cores
        pool = Pool(n_cores)

        #set up async processes
        processes = [pool.apply_async(coconut, args=(index, x))
                     for x in range(50)]
        
        #run processes
        results = [proc.get() for proc in processes]
        


