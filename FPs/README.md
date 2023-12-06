# Folder structure

- **coconut**: Is a placeholder for containing the fingerprints calculated by `fp_script.py` for the **COCONUT**: database.  

- **drug_rep_hub**: Is a placeholder for containing the fingerprints calculated by `fp_script.py` for the **drug_rep_hub** database.  

- **antimicrobial**: Is a placeholder for containing the fingerprints calculated by `make_clf_script.py` for the antimicrobial activity prediction dataset.  

- **antiviral**: Is a placeholder for containing the fingerprints calculated by `make_clf_script.py` for the antiviral activity prediction dataset.  

- **cytotoxic**: Is a placeholder for containing the fingerprints calculated by `make_clf_script.py` for the cytotoxic activity prediction dataset.  

Each fingerprint is saved as e.g. "ecfp.pkl" in the respective folder. The files themselves cannot be uploaded due to space limitations (the one inside is simply a placeholder), but the calculation can be reproduced by processing `coconut.csv` or `drug_rep_hub.csv` with `fp_script.py`, or running `make_clf_script.py` with a given property (e.g. "antimicrobial").  

