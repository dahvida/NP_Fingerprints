"""Classification analysis script.

Script to train Random Forest (RF) models for each taxonomy class from
the COCONUT DB, using all computed fingerprints.

Before running this script, make sure there are .pkl files in the
../FPs folder. If not, run the fp_script.py file.

Each Random Forest model is optimized on the validation set and evaluated
on the test set. The training-validation-test splits are done via 
scaffold splitting with an 80:10:10 ratio. Model performance is expressed
as:
    - Balanced accuracy
    - Matthews correlation coefficient
    - F1 score
    - Precision
    - Recall
    - ROC-AUC
    - PR-AUC

The procedure goes as follows:
    1. Load labels and create scaffold split indices. Indices will be
       saved for reproducibility in ../Results
    2. Check that the splits don't have excessive class imbalance
    3. Loop over each of the 5 taxonomies available:
        3.1. Loop over each precomputed fingerprint:
            3.1.1. Load and split fingerprints&labels in train-val-test
            3.1.2. Optimize RF model (if classifier == "RF")
            3.1.3. Train and evaluate on test set using optimal params
            3.1.4. Store evaluation results
        3.2. Save results with all fingerprints for a given class as
             ../Results/classification/class_name.csv
    
Users can tune the number of optimization iterations, the classification
metric to optimize, the number of replicates for performance evaluation
and the random seed for scaffold splitting.
"""

from clf import *
from utils import *
import pandas as pd
import numpy as np
import argparse
import os

######################################################################

parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--opt_iterations",
                    default = 20,
                    type = int,
                    help = "number of iterations during optimization")

parser.add_argument("--opt_metric",
                    default = "PR-AUC",
                    help = "metric to optimize during optimization")

parser.add_argument("--n_replicates",
                    default = 5,
                    type = int,
                    help = "number of replicates for performance evaluation")

parser.add_argument("--classifier",
                    default = "RF",
                    help = "type of classifier [RF, NB]")

parser.add_argument('--random_seed',
                    default = 42,
                    type = int,
                    help = "random seed for reproducibility")

args = parser.parse_args()

######################################################################

def main(
        opt_iterations,
        opt_metric,
        n_replicates,
        classifier,
        random_seed
        ):
    
    #print run params
    print("[clf]: Starting run with parameters:")
    print(f"-- optimization iterations: {opt_iterations}")
    print(f"-- optimization metric: {opt_metric}")
    print(f"-- evaluation replicates: {n_replicates}")
    print(f"-- classifier: {classifier}")
    print(f"-- random seed: {random_seed}")
    
    #create scaffold split and save in ../Results
    print("[clf]: Generating scaffold split")
    df = pd.read_csv("../Data/clean_coconut.csv")
    train_idx, val_idx, test_idx = split(df, random_seed)
    pickle_save(train_idx, "../Results/train_idx.pkl", verbose=False)
    pickle_save(val_idx, "../Results/val_idx.pkl", verbose=False)
    pickle_save(test_idx, "../Results/test_idx.pkl", verbose=False)
    print("[clf]: Splits generated, indexes saved in ../Results")
    
    #Check that the splits are not affected by too much class imbalance.
    #especially relevant for i.e. marine class (0.8% of dataset)
    print("[clf]: Loading labels")
    labels = pickle_load("../Data/labels.pkl")
    assert_splits(labels, train_idx, val_idx, test_idx)
    print("[clf]: All splits passed quality check for number of actives")

    #get FP names and paths to pickle them
    names = os.listdir("../FPs")
    fp_names = [x[:-4] for x in names]
    fp_paths = ["../FPs/" + x for x in names] 
    class_names = ["bacteria", "plant", "fungi", "animal", "marine"]
    
    #loop over each class
    for i in range(labels.shape[1]):
        print("-----------------------------------------")
        print(f"[clf]: Analysing class {class_names[i]}")
        
        #slice labels for the i-th class
        label_idx = labels[:,i]
        y_train = label_idx[train_idx]
        y_val = label_idx[val_idx]
        y_test = label_idx[test_idx]

        #create results container
        results_list = []
        
        #loop over all available fingerprints
        for j in range(len(fp_paths)):
            #load j-th FP and split it in train-val-test
            fp = pickle_load(fp_paths[j], verbose=False)
            X_train = fp[train_idx,:]
            X_val = fp[val_idx, :]
            X_test = fp[test_idx, :]
            
            if classifier == "RF":
                #optimize RF
                print(f"[clf]: Starting optimization for {fp_names[j]}")
                optimum = optimize(X_train,
                       y_train,
                       X_val,
                       y_val,
                       opt_iterations,
                       opt_metric,
                       )
            else:
                optimum = None

            results = evaluate(X_train,
                       y_train,
                       X_test,
                       y_test,
                       optimum,
                       classifier,
                       n_replicates)
            

            #store results in container
            results_list.append(results)
        
        #save results for all FPs for the i-th class
        store_results(
                results_list,
                fp_names,
                class_names[i],
                algorithm
                )
        print("[clf]: Results saved at ../Results/classification")

    

if __name__ == "__main__":
    main(
        args.opt_iterations,
        args.opt_metric,
        args.n_replicates,
        args.classifier,
        args.random_seed
        )






