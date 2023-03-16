# Folder structure

- [classification:](clean_coconut.csv) Contains classification performance for each calculated fingerprint. Results are summarized in terms of balanced accuracy, matthews correlation coefficient, F1 score, precision, recall, ROC-AUC and PR-AUC.  

- [correlation_stats:](correlation_stats.pkl) Contains the correlation matrix between
fingerprints for each sampling iteration.  

- [sample_stats:](sample_stats.txt) Contains the dataset statistics for each batch sampled
from the preprocessed COCONUT database.  

- [similarity_stats:](similarity_stats.csv) Contains distribution statistics for the similarities during each sampling iteration according to each fingerprint.  

- [index.pkl:](index.pkl) Matrix defining which compounds to sample during each iteration.  

- [test_idx.pkl:](test_idx.pkl) Index of the compounds used as test set.  

- [train_idx.pkl:](train_idx.pkl) Index of compounds used as training set.  

- [val_idx.pkl:](val_idx.pkl) Index of compounds used as validation set.  

