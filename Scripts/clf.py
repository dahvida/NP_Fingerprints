"""
Utility functions used in clf_script.py
"""

from sklearn.metrics import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import *
from sklearn.preprocessing import StandardScaler
from hyperopt import tpe, hp, fmin, Trials
import numpy as np
import deepchem as dc
import pandas as pd
from typing import *

######################################################################

def optimize(
        x_train: np.ndarray,
        y_train: np.ndarray,
        x_val: np.ndarray,
        y_val: np.ndarray,
        iters: int = 20,
        metric: str = "PR-AUC"
            ) -> Dict:
    """Hyperparameter optimization function
    
    Find optimal values for number of estimators, maximum tree depth,
    minimum samples per split, minimum samples per leaf and number of
    features per split according to the chosen optimization metric.

    Args:
        x_train:    (M, K) training set fingerprints
        y_train:    (M, 1) training set labels
        x_val:      (N, K) validation set fingerprints
        y_val:      (N, 1) validation set fingerprints
        iters:      number of iterations for optimization
        metric:     performance metric to optimize (either PR-AUC
                    or ROC-AUC)

    Returns:
        A dictionary containing the optimal hyperparameter
        values
    """        
            
    #load correct metric depending on input
    if metric == "PR-AUC":
        metric = average_precision_score
    else:
        metric = roc_auc_score
            
    #define optimization grid
    space = {
        'n_estimators':         hp.quniform('n_estimators',
                                            50, 500, 50),

        'max_depth':            hp.quniform('max_depth',
                                            5, 15, 2),

        'min_samples_split':    hp.quniform('min_samples_split',
                                            2, 20, 1),

        'min_samples_leaf':     hp.quniform('min_samples_leaf',
                                            2, 100, 1),

        'max_features':         hp.choice('max_features',
                                            ["sqrt", "log2", 0.1])
            }

    #define eval function
    def model_eval(args):
                
        #parse vars into int for RF model
        args["n_estimators"] = int(args["n_estimators"])
        args["max_depth"] = int(args["max_depth"])
        args["min_samples_split"] = int(args["min_samples_split"])
        args["min_samples_leaf"] = int(args["min_samples_leaf"])

        #define model with custom args, using all cores
        model = RandomForestClassifier(
                        class_weight = "balanced_subsample",
                        n_jobs = -1,
                        **args
                        )
                
        #fit model and return logits
        model.fit(x_train, y_train)
        probs = model.predict_proba(x_val)[:,1]
                                
        return 1 - metric(y_val, probs)
            
    #create trials object
    trials = Trials()
            
    #run optimization
    optimum = fmin(
                fn = model_eval,
                space = space,
                algo = tpe.suggest,
                max_evals = iters,    
                trials = trials,
                verbose = False
                )
            
    #parse optimum vars into ints or into correct value
    optimum['max_features'] = ["sqrt", "log2", 0.1][optimum['max_features']]
    optimum["n_estimators"] = int(optimum["n_estimators"])
    optimum["max_depth"] = int(optimum["max_depth"])
    optimum["min_samples_split"] = int(optimum["min_samples_split"])
    optimum["min_samples_leaf"] = int(optimum["min_samples_leaf"]) 

    return optimum

#--------------------------------------------------------------------------#

def evaluate(
        x_train: np.ndarray,
        y_train: np.ndarray,
        x_test: np.ndarray,
        y_test: np.ndarray,
        params: Dict,
        algorithm: str,
        iters: int = 5
             ) -> np.ndarray:
    """Model performance evaluation function

    Args:
        x_train:    (M, K) training set fingerprints
        y_train:    (M, 1) training set labels
        x_test:     (J, K) test set fingerprints
        y_test:     (J, 1) test set labels
        params:     parameters to use to construct the RF model
        algorithm:  toggles whether to use RF or NB as classifier
        iters:      number of training and evaluation repetitions
    
    Returns:
        A numpy array (7,2) containing in the first column the mean 
        performance for each metric across iterations, in the second column
        their standard deviation.
        The metrics are ordered as follows:
            - Balanced accuracy
            - Matthews correlation coefficient
            - F1 Score
            - Precision
            - Recall
            - ROC-AUC
            - PR-AUC
    """
    
    #create containers with right dimensions
    storage = np.zeros((7, iters))
    results = np.zeros((7, 2))
    
    #loop analysis over number of iterations
    for i in range(iters):
        
        if algorithm == "RF":
        #create model with param config
            model = RandomForestClassifier(
                class_weight = "balanced_subsample",
                n_jobs = -1,
                **params
                )
        else:
            model = GaussianNB()
            scaler = StandardScaler()
            x_train = scaler.fit_transform(x_train)
            x_test = scaler.transform(x_test)

        #fit model on training data
        model.fit(x_train, y_train)

        #get binary predictions and logits for test set
        pred_label = model.predict(x_test)
        pred_prob = model.predict_proba(x_test)[:,1]
        
        #measure all metrics and store them in storage array
        storage[0,i] = balanced_accuracy_score(y_test, pred_label)
        storage[1,i] = matthews_corrcoef(y_test, pred_label)
        storage[2,i] = f1_score(y_test, pred_label)
        storage[3,i] = precision_score(y_test, pred_label)
        storage[4,i] = recall_score(y_test, pred_label)
        storage[5,i] = roc_auc_score(y_test, pred_prob)
        storage[6,i] = average_precision_score(y_test, pred_prob)
    
    #get mean and std across repetition and store them in output array
    results[:,0] = np.mean(storage, axis=1)
    results[:,1] = np.std(storage, axis=1)

    return results

#--------------------------------------------------------------------------#

def store_results(
        results: List[np.ndarray],
        fp_names: List[str],
        class_name: str,
        algorithm: str
        ) -> None:
    """Converts list of results in .csv file
    
    Args:
        results:    List (X) of outputs from evaluate function, each entry
                    represents the results from a given fingerprint
        fp_names:   List (X) of names of fingerprints used
        class_name: name of class used for the evaluation analysis
        algorithm:  name of the classification algorithm used for the analysis

    Returns:
        None
    """

    #create container for all results and indexes for correct allocation
    array = np.zeros((len(results), 14))
    idx_means = np.array([0,2,4,6,8,10,12])
    idx_std = np.array([1,3,5,7,9,11,13])
    
    #store results in container as mean1, std1, mean2, std2, mean3, std3 ...
    for i in range(len(results)):
        array[i, idx_means] = results[i][:,0]
        array[i, idx_std] = results[i][:,1]
    
    #create metric name template
    metric_names = [
            "Balanced accuracy",
            "MCC",
            "F1",
            "Precision",
            "Recall",
            "ROC-AUC",
            "PR-AUC"
            ]
    
    #create proper names by modifying common template
    means = [x + " - mean" for x in metric_names]
    std = [x + " - std" for x in metric_names]

    #interleave lists so that the output is mean1, std1, mean2, std2 ...
    metric_names = [val for pair in zip(means, std) for val in pair]
    
    #create dataframe and save in correct directory
    df = pd.DataFrame(
            data = array,
            index = fp_names,
            columns = metric_names
            )
    
    if algorithm == "RF:
    	df.to_csv("../Results/classification_rf/" + class_name + ".csv")
    else:
    	df.to_csv("../Results/classification_nb/" + class_name + ".csv")

#--------------------------------------------------------------------------#

def split(
        df: pd.DataFrame,
        seed: int = 0
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scaffold splitting function
    
    Args:
        df:     dataset to split, must have a column named "sugar_free_smiles"
        seed:   random seed to make splits reproducible

    Returns
        Tuple of 3 numpy arrays containing the sample indexes to split
        into training, validation and test sets. Split is done with 80:10:10
        ratio.
    """

    #fetch smiles from df as list
    smiles = list(df["sugar_free_smiles"])
    
    #create splitter object
    splitter = dc.splits.ScaffoldSplitter()
    
    #create empty deepchem dataset object with smiles
    dummy = np.zeros(len(smiles))
    dataset = dc.data.DiskDataset.from_numpy(X = dummy,
                                             y = dummy,
                                             w = dummy,
                                             ids = smiles)
    
    #split deepchem dataset
    train_idx, val_idx, test_idx = splitter.split(dataset, seed=seed)
    
    return np.array(train_idx), np.array(val_idx), np.array(test_idx)

def assert_splits(
        labels: np.ndarray,
        train_idx: np.ndarray,
        val_idx: np.ndarray,
        test_idx: np.ndarray
        ) -> None:
    """Quality check for class distribution in splits
    
    Since scaffold split is employed, some classes might not have a 80:10:10
    distribution between train-val-test sets. The function ensures that the
    chosen random seed doesn't cause excessive class imbalance in the sets.
    
    Args:
        labels:         (A, 5) array with class labels for all compounds
        train_idx:      (M, 1) array with IDs of training set componds
        val_idx:        (N, 1) array with IDs of validation set compounds
        test_idx:       (J, 1) array with IDs of test set compounds

    Returns:
        None
    """
    
    for i in range(labels.shape[1]):
        y_train = labels[train_idx, i]
        y_val = labels[val_idx, i]
        y_test = labels[test_idx, i]
        assert np.sum(y_train) > 500, "[clf]: Split failed quality check, try different random seed"
        assert np.sum(y_val) > 50,  "[clf]: Split failed quality check, try different random seed"
        assert np.sum(y_test) > 50, "[clf]: Split failed quality check, try different random seed"











