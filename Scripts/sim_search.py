"""
Utiliy functions for sim_search_script.py
"""

from scipy.spatial.distance import jaccard
import numpy as np
import pandas as pd
from typing import *

#############################################################################

def sim_search(
        array: np.ndarray,
        sim_metric: str = "default"
        ) -> np.ndarray:
    """Runs efficient similarity search on a given array
    
    Args:
        array:      matrix (M,K) to use for similarity search
        sim_metric: metric type to use for the search

    Returns:
        A vector ((M**2 - M) / 2) of all unique pairwise similarities
        for each element in the array
    """
    
    fp = array.copy()

    #swap metric type (used for MHFP & MAP4)
    if sim_metric == "default":
        metric = jaccard
        #fix array for count-based fingerprints (AP, TT, Avalon)
        fp[fp > 1] = 1
    else:
        metric = jaccard_like
    
    #preallocate vector of correct size
    n_samples = fp.shape[0]
    n_comps = int((n_samples ** 2 - n_samples) / 2)
    similarity = np.zeros((1, n_comps,), dtype=np.float32)
    
    #start pairwise sim counter
    count = 0
    
    #loop over all unique pairs (skip redundant comparisons)
    for i in range(n_samples - 1):
        for j in range(i+1, n_samples):
            #calculate sim as 1 - jaccard distance
            similarity[0, count] = 1 - metric(fp[i,:],
                                            fp[j,:])
            count += 1

    return similarity

#--------------------------------------------------------------------------#

def jaccard_like(a,b):
    """
    Used instead of normal tanimoto for MAP4 and MHFP
    """
    return 1 - float(np.count_nonzero(a == b)) / float(len(a))

#--------------------------------------------------------------------------#

def eval_sim(
        similarity: np.ndarray
        ) -> np.ndarray:
    """Collect all statistics for a similarity vector
    
    Args:
        similarity:     similarity vector (S,) to analyse

    Returns:
        A vector (22,) containing all statistical indexes describing the original
        similarity vector
    """
    
    #preallocate
    output = np.zeros((22))
    
    #get mean, standard deviation and median
    output[0] = np.mean(similarity)
    output[1] = np.std(similarity)
    output[2] = np.median(similarity)
    
    #get percentiles with an interval of 5 (i.e. 5,10,15,20,25...)
    percentiles = np.arange(5, 100, 5)
    output[3:] = np.percentile(similarity, percentiles)
    
    return output

#--------------------------------------------------------------------------#

def save_sim_df(
        output_box: np.ndarray,
        fp_names: List[str],
        path: str,
        verbose: bool = True
        ) -> None:
    """Saves similarity stats array as csv
    
    Args:
        output_box:     matrix (S,20) of all similarity vectors
        fp_names:       fingerprint names
        path:           path where to save the .csv file
        verbose:        whether to print the saving directory

    Returns:
        None
    """

    #create names for all collected statistics
    index_names = []
    index_names.append("Mean")
    index_names.append("STD")
    index_names.append("Median")
    percentiles = np.arange(5, 100, 5)
    for i in range(len(percentiles)):
        index_names.append("Percentile: " +
                           str(percentiles[i]))
    
    #create and save dataframe
    df = pd.DataFrame(
            data = output_box,
            index = index_names,
            columns = fp_names
            )
    df.to_csv(path)
    if verbose is True:
        print(f"[sim_search]: Saving similarity stats as {path}")

#--------------------------------------------------------------------------#

def save_corr_df(
        correlation_matrix: np.ndarray,
        fp_names: List[str],
        path: str,
        verbose: bool = True
        ) -> None:
    """Saves correlation matrix as csv
    
    Args:
        correlation_matrix:     correlation matrix (20,20) between all
                                fingerprints
        fp_names:               fingerprint names
        path:                   path where to save the .csv file
        verbose:                whether to print the saving directory

    Returns:
        None
    """

    #create and save dataframe
    df = pd.DataFrame(
            data = correlation_matrix,
            index = fp_names,
            columns = fp_names
            )
    df.to_csv(path)
    if verbose is True:
        print(f"[sim_search]: Saving correlation matrix as {path}")






