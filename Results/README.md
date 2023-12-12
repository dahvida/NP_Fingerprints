# Folder structure

- **antibiotic:** Contains scaffold split indices and classification performancefor each calculated fingerprint using RF and MLP on the antibiotic dataset. Results are summarized in terms of matthews correlation coefficient, ROC-AUC, PR-AUC, precision, recall and specificity.  

- **antiviral:** Same as above but for the antiviral dataset.  

- **cytotoxic:** Same as above but for the cytotoxicity dataset.  

- **antileishmanial:** Same as above but for the antileishmanial dataset.  

- **antifungal:** Same as above but for the antifungal dataset.  

- **antitumour:** Same as above but for the antitumour dataset.  

- **hiv:** Same as above but for the hiv dataset.  

- **kinase_c:** Same as above but for the kinase_c dataset.  

- **serine_protease:** Same as above but for the serine_protease dataset.  

- **anti_inflammatory:** Same as above but for the anti_inflammatory dataset.  

- **phosphatase:** Same as above but for the phosphatase dataset.  

- **atpase:** Same as above but for the atpase dataset.  

- **antimalarial:** Same as above but for the antimalarial dataset.  

- **correlation_stats:** Contains the correlation matrix between fingerprints for each sampling iteration.  

- **sample_stats:** Contains the dataset statistics for each batch sampled from the preprocessed COCONUT database.  

- **similarity_stats:** Contains distribution statistics for the similarities during each sampling iteration according to each fingerprint.  

- **drug_rep_hub:** Contains the correlation matrix and similarity distribution for the Drug Repurposing Hub dataset.  

- **index.pkl:** Matrix defining which compounds to sample during each iteration.  

