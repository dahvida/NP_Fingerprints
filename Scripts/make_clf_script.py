"""Classification dataset generation script.

Script to generate a classification dataset from the labels provided by the
CMNPD database of marine natural products. The script generates all folders
and files necessary to then execute clf_script.py to evaluate the performance
of all fingerprints on the dataset.

The procedure goes as follows:
    1. Select all molecules with the desired property, i.e. "antimicrobial"
    2. For the negative class, select either an equal number of random NPs
       from CMNPD, or enough so that the dataset size is equal to "min_dataset_size"
    3. Save the generated dataset
    4. Compute all fingerprints for the dataset

Users can choose the property to select, how to name the generated dataset
and the minimum size of the dataset, in case there is only a handful of samples
with the desired property.
"""

import pandas as pd
import argparse
import os

###########################################################################

parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--prop_name",
                    help = "Property to extract from labeled NPs from CMNPD")

parser.add_argument("--output_name",
                    help = "Name to use when saving the generated dataset")

parser.add_argument("--min_dataset_size",
                    default = 1000,
                    type = int,
                    help = "Minimum size of the dataset")

args = parser.parse_args()

###########################################################################

def main(prop_name,
         output_name,
         min_dataset_size):
    
    print("[make_clf]: Starting dataset generation process")
    
    # make all necessary folders where to save the data
    fp_path = "../FPs/" + output_name
    clf_path = "../Results/" + output_name
    db_path = "../Data/" + output_name + ".csv"
    if not os.path.exists(fp_path):
        os.makedirs(fp_path)
        os.makedirs(clf_path)
    
    # load dataframe
    db = pd.read_csv("../Data/cmnpd.csv", index_col=0)
    
    # get relevant slices
    db_prop = db[db['BIOACTIVITY_IN_BRIEF'].str.contains(prop_name, case=False, na=False)]
    db_no_prop = db[~db['BIOACTIVITY_IN_BRIEF'].str.contains(prop_name, case=False, na=False)]
    db_no_prop.sample(frac=1)
    
    # get either an equal number of inactives (randomly) or a number so that
    # n_act + n_inact = min_dataset_size
    if len(db_prop) < min_dataset_size/2:
        t = min_dataset_size - len(db_prop)
    else:
        t = len(db_prop)
    
    db_no_prop = db_no_prop.head(t)
    
    labels = [1]*len(db_prop) + [0]*len(db_no_prop)
    
    print(f"[make_clf]: Dataset generated with {len(db_prop)} actives and {len(db_no_prop)} sampled inactives")
    
    db_final = pd.concat([db_prop, db_no_prop], axis=0)
    db_final["y"] = labels
    db_final = db_final.sample(frac=1)
    
    db_final.to_csv(db_path)
    
    print("[make_clf]: Starting fingerprint generation...")
    os.system("python3 fp_script.py --FP_type rdkit --dataset " + output_name )
    os.system("python3 fp_script.py --FP_type jmap --dataset " + output_name )
    os.system("python3 fp_script.py --FP_type cdk --dataset " + output_name )
    os.system("python3 fp_script.py --FP_type minhash --dataset " + output_name )
    
if __name__ == "__main__":
    main(
        args.prop_name,
        args.output_name,
        args.min_dataset_size
    )

    
    
    
