"""Similarity search script.

Script to calculate pairwise similarities between natural products
from the COCONUT DB.

Before running this script, make sure there are .pkl files in the
../FPs folder. If not, run the fp_script.py file.

Due to dataset size, it is impossible to compute all pairwise
similarities in one go. As such, the analysis is repeated across
randomly sampled batches.

The procedure goes as follows:
    1. Generate indices to use to sample from COCONUT DB
    2. For each iteration:
        2.1. Load dataset sample and collect statistics (class
             frequencies and scaffold diversity)
        2.2. Load i-th fingerprint, slice the array according to
             indices and compute similarities
        2.3. Repeat step 2.2 for all fingerprints
        2.4. Store correlations and similarity statistics
    3. Save correlations to ../Results/correlation_stats, 
       similarity statistics to ../Results/similarity_stats,
       sample statistics to ../Results/sample_stats and indices
       for sampling as ../Results/index.pkl

Parallelization is implemented so that each batch iteration can
be run in parallel on a different core. Keep in mind that each process
can use several GBs of RAM, as such the limiting factor is likely
to be the total RAM rather than the number of cores.

Users can tune how many compounds to sample, how many sampling iterations
to do, how many cores to use and which random seed to set for generating
sampling indices
"""

from sim_search import *
from utils import *
import pandas as pd
import numpy as np
import os
from multiprocessing import Process, Pool
import argparse

######################################################################

parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--sample_size",
                    default = 10000,
                    type = int,
                    help = "number of samples to draw for each batch")

parser.add_argument("--n_replicates",
                    default = 50,
                    type = int,
                    help = "number of replicates to use for the analysis")

parser.add_argument("--n_cores",
                    default = 20,
                    type = int,
                    help = "number of cores to use for the analysis")

parser.add_argument('--random_seed',
                    default = 42,
                    type = int,
                    help = "random seed for reproducibility")

args = parser.parse_args()

######################################################################

def main(index, iteration_n):
    
    #get relevant index according to iteration number
    print(f"[sim_search]: Running iteration {iteration_n}")
    index = index[:,iteration_n]
    n_compounds = index.shape[0]
    
    #get names of all FPs available in ../FPs
    names = os.listdir("../FPs/coconut")
    fp_names = [x[:-4] for x in names]
    fp_paths = ["../FPs/coconut" + x for x in names]
    
    #compute number of pairwise similarity calculations and
    #preallocate array of correct size
    total_pairs = int((n_compounds**2 - n_compounds) / 2 )
    pair_matrix = np.zeros((len(fp_paths), total_pairs),
                           dtype = np.float32)
    stats_matrix = np.zeros((22, len(fp_paths)),
                             dtype = np.float32)
    
    #create paths to save files
    path_sim = "../Results/similarity_stats/" 
    path_corr = "../Results/correlation_stats/"
    path_df = "../Results/sample_stats/"
    
    #load dataset slice i-th and collect stats
    db = pd.read_csv("../Data/coconut.csv")
    db = db.iloc[index]
    db_stats = get_statistics(db)
    db_stats.to_csv(path_df + str(iteration_n) + ".csv")
    
    #loop over each FP in ../FPs
    for j in range(len(fp_paths)):
        
        #load FP
        fp = pickle_load(fp_paths[j], verbose=False)
        
        #overwrite FP with its slice according to index (important to 
        #limit memory usage to the bare minimum we need for i-th run)
        fp = fp[index,:]

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
    correlation_matrix = np.corrcoef(pair_matrix)
    
    #save everything
    save_corr_df(correlation_matrix, fp_names,
                     path = path_corr + str(iteration_n) + ".csv",
                     verbose = False)
    save_sim_df(stats_matrix, fp_names,
                    path = path_sim + str(iteration_n) + ".csv",
                    verbose = False)
    


if __name__ == "__main__":
    
    #load user arguments from parser
    seed = args.random_seed
    sample_size = args.sample_size
    n_replicates = args.n_replicates
    n_cores = args.n_cores

    #set seed for index generation
    print(f"[sim_search]: Beginning analysis")
    np.random.seed(seed)

    #generate matrix to sample from the cleaned dataset for the
    #correlation analysis
    index = np.random.randint(low = 0,
                              high = 129869,
                              size = (sample_size, n_replicates))

    #save as pkl in ../Results
    pickle_save(index, "../Results/index.pkl", verbose = False)

    #write run info
    n_pairs = int((sample_size ** 2 - sample_size) / 2)
    print(f"[sim_search]: Run settings:")
    print(f"-- compounds: {sample_size}")
    print(f"-- similarity pairs: {n_pairs}")
    print(f"-- iterations: {n_replicates}")
    print(f"-- cores: {n_cores}")
    print(f"-- random seed: {seed}")
    print(f"-- Using smooth tanimoto for AP, TT, Avalon")
    
    #create pool with n cores
    pool = Pool(n_cores)

    #set up async processes
    processes = [pool.apply_async(main, args=(index, x)) for x in range(n_replicates)]
    
    #run processes
    results = [proc.get() for proc in processes]

    print(f"[sim_search]: Similarity statistics have been saved in ../Results/similarity_stats")
    print(f"[sim_search]: Correlation matrix has been saved in ../Results/correlation_stats")
    print(f"[sim_search]: Sample statistics have been saved in ../Results/sample_stats")


