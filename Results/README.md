# Folder structure

- **antimicrobial:** Contains scaffold split indices and classification performancefor each calculated fingerprint using RF and MLP on the antimicrobial dataset. Results are summarized in terms of matthews correlation coefficient, ROC-AUC, PR-AUC, precision, recall and specificity.  

- **antiviral:** Same as above but for the antiviral prediction dataset.  

- **cytotoxic:** Same as above but for the cytotoxicity prediction dataset.  

- **correlation_stats:** Contains the correlation matrix between fingerprints for each sampling iteration.  

- **sample_stats:** Contains the dataset statistics for each batch sampled from the preprocessed COCONUT database.  

- **similarity_stats:** Contains distribution statistics for the similarities during each sampling iteration according to each fingerprint.  

- **drug_rep_hub:** Contains the correlation matrix and similarity distribution for the Drug Repurposing Hub dataset.  

- **index.pkl:** Matrix defining which compounds to sample during each iteration.  

