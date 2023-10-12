# Folder structure

- **classification:** Contains classification performance for each calculated fingerprint. Results are summarized in terms of balanced accuracy, matthews correlation coefficient, F1 score, precision, recall, ROC-AUC and PR-AUC.  

- **correlation_stats:** Contains the correlation matrix between
fingerprints for each sampling iteration.  

- **sample_stats:** Contains the dataset statistics for each batch sampled
from the preprocessed COCONUT database.  

- **similarity_stats:** Contains distribution statistics for the similarities during each sampling iteration according to each fingerprint.  

- **drug_rep_hub:** Contains the correlation matrix and similarity distribution for the Drug Repurposing Hub dataset.  

- **index.pkl:** Matrix defining which compounds to sample during each iteration.  

- **test_idx.pkl:** Index of the compounds used as test set.  

- **train_idx.pkl:** Index of compounds used as training set.  

- **val_idx.pkl:** Index of compounds used as validation set.  

