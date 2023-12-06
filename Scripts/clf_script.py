"""Classification analysis script.

Script to train Random Forest (RF) or Multilayer Perceptron (MLP) models 
for a classification task from the CMNPD database.

Before running this script, generate the desired classification task using
make_clf_script.py.

Each model is optimized on the validation set and evaluated
on the test set. The training-validation-test splits are done via 
scaffold splitting with an 80:10:10 ratio. Model performance is expressed
as:
    - Matthews correlation coefficient
    - Precision
    - Recall
    - Specificity
    - ROC-AUC
    - PR-AUC

The procedure goes as follows:
    1. Load labels and create scaffold split indices. Indices will be
       saved for reproducibility in ../Results
    2. Check that the splits don't have excessive class imbalance
    3. Loop over each precomputed fingerprint:
            3.1. Load and split fingerprints&labels in train-val-test
            3.2. Optimize model 
            3.3. Train and evaluate on test set using optimal params
            3.4. Store evaluation results
    4. Save results with all fingerprints for a given class as
             ../Results/dataset_name/classification_algorithm.csv
    
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
                    default = "ROC-AUC",
                    help = "metric to optimize during optimization")

parser.add_argument("--n_replicates",
                    default = 5,
                    type = int,
                    help = "number of replicates for performance evaluation")

parser.add_argument("--classifier",
                    default = "RF",
                    help = "type of classifier [RF, MLP]")

parser.add_argument("--dataset",
                    help = "dataset to process [antibiotic_np, drug_vs_np]")

args = parser.parse_args()

######################################################################

def main(
        opt_iterations,
        opt_metric,
        n_replicates,
        classifier,
        dataset
        ):
    
    #print run params
    print("[clf]: Starting run with parameters:")
    print(f"-- optimization iterations: {opt_iterations}")
    print(f"-- optimization metric: {opt_metric}")
    print(f"-- evaluation replicates: {n_replicates}")
    print(f"-- classifier: {classifier}")
    
    # load data
    df = pd.read_csv("../Data/" + dataset + ".csv")
    
    #create scaffold split and save in ../Results/dataset, otherwise skip
    if not os.path.exists("../Results/" + dataset + "/train_idx.pkl"):
        print("[clf]: Generating scaffold split")
        train_idx, val_idx, test_idx = split(df)
        pickle_save(train_idx, "../Results/" + dataset + "/train_idx.pkl", verbose=False)
        pickle_save(val_idx, "../Results/" + dataset + "/val_idx.pkl", verbose=False)
        pickle_save(test_idx, "../Results/" + dataset + "/test_idx.pkl", verbose=False)
        print(f"[clf]: Splits generated, indexes saved in ../Results/{dataset}")
    else:
        print(f"[clf]: Using precomputed splits found in ../Results/{dataset}")
        train_idx = pickle_load("../Results/" + dataset + "/train_idx.pkl", verbose=False)
        val_idx = pickle_load("../Results/" + dataset + "/val_idx.pkl", verbose=False)
        test_idx = pickle_load("../Results/" + dataset + "/test_idx.pkl", verbose=False)
    
    #Check that the splits are not affected by too much class imbalance.
    labels = np.array(df["y"])
    assert_splits(labels, train_idx, val_idx, test_idx)
    print("[clf]: All splits passed quality check for number of actives")

    #get FP names and paths to pickle them
    names = os.listdir("../FPs/" + dataset + "/")
    fp_names = [x[:-4] for x in names]
    fp_paths = ["../FPs/" + dataset + "/" + x for x in names] 
    
    #loop over each class
    print("-----------------------------------------")
    print(f"[clf]: Starting classification evaluation")
        
    #slice labels for the i-th class
    y_train = labels[train_idx]
    y_val = labels[val_idx]
    y_test = labels[test_idx]

    #create results container
    results_list = []
        
    #loop over all available fingerprints
    for j in range(len(fp_paths)):
            #load j-th FP and split it in train-val-test
            fp = pickle_load(fp_paths[j], verbose=False)
            X_train = fp[train_idx,:]
            X_val = fp[val_idx, :]
            X_test = fp[test_idx, :]
            
            #optimize hyperparams
            print(f"[clf]: Starting optimization for {fp_names[j]}")
            optimum = optimize(X_train,
                       y_train,
                       X_val,
                       y_val,
                       opt_iterations,
                       opt_metric,
                       classifier
                       )

            #eval with n_replicates
            print(f"[clf]: Starting evaluation for {fp_names[j]}")
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
                classifier,
                dataset
                )
    print(f"[clf]: Results saved at ../Results/{dataset}")

    

if __name__ == "__main__":
    main(
        args.opt_iterations,
        args.opt_metric,
        args.n_replicates,
        args.classifier,
        args.dataset
        )






