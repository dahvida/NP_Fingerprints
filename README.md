# NP_Fingerprints
![Python 3.7](https://img.shields.io/badge/python-3.7%20%7C%203.8-brightgreen)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)  
Scripts to calculate fingerprints and simiilarity matrices for natural product databases.  

## Repository structure
- [Data:](Data) Contains the COCONUT database, its preprocessed version and dataset statistics.  
- [FPs:](FPs) Contains all calculated fingerprints as numpy arrays, stored via pickle.  
- [Scripts:](Scripts) Contains all scripts and utility functions used to calculate fingerprints and reproduce the results.  
- [Results:](Results) Contains all the results from the analysis.  

## Installation  
All necessary Python packages can be installed via conda from the `environment.yml` file.  
```
git clone https://github.com/dahvida/NP_Fingerprints
conda env create --name np_fp --file=environment.yml
conda activate np_fp
cd NP_Fingerprints/Scripts
```
Additionally, you need to install Java 11 for computing CDK / jCompoundMapper fingerprints, and unzip the .jar files in [FP_calc](Scripts/FP_calc).  

## Tutorial
**All scripts must be executed from the [Scripts](Scripts) folder.**  
**Make sure to unzip `COCONUT_DB.zip` in [Data](Data) and the Java tools in [Scripts/FP_calc](Scripts/FP_calc).**  
To clean the raw COCONUT database (accessed from https://github.com/reymond-group/Coconut-TMAP-SVM), use `cleanup_script.py`. The
preprocessed dataset will be saved in [Data](Data), as well as class and Murcko scaffold statistics:  
```
python3 cleanup_script.py
```
To calculate a set of fingerprints for the preprocessed version of the COCONUT database, use `fp_script.py`. You must specify which
set of fingerprints you're interested on calculating via the `--FP_type` flag. Options are "rdkit", "minhash", "cdk" and "jmap". All
fingerprints will be stored in pickle files as numpy arrays in the [FPs](FPs) folder.  
```
python3 fp_script.py --FP_type rdkit
```
To generate the fingerprint correlation matrices, use `sim_search_script.py`. The script will sample 50 batches of 10.000 compounds each
to reduce the computational load of running all pairwise comparisons on 129k compounds. The output will be saved in the [Results](Results)
folder. The default arguments are the same as the ones
employed in our study. Otherwise, they can be modified via the appropriate flags, i.e:  
```
python3 sim_search_script.py --sample_size 5000 --n_cores 4
```
To generate the classification analysis results, use `clf_script.py`. By default, the script will run 20 bayesian hyperparameter optimization
iterations and evaluate its performance on the test set for 5 replicates. Similarly as before, these parameters can be modified via the 
appropriate flags, i.e:   
```
python3 clf_script.py --opt_iters 10 --n_replicates 50
```
Finally, to get more information on the arguments that can be passed to each command line tool and the steps employed in each script, you can
use the `--help` flag:
```
python3 clf_script.py --help
```

## Further information
More details on which fingerprints are available, which results are saved where and so forth can be found in the respective README files of
each folder.  

## How to cite
Link to publication goes here


