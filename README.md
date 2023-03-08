# NP_Fingerprints
![Python 3.7](https://img.shields.io/badge/python-3.7%20%7C%203.8-brightgreen)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)  
Scripts to calculate fingerprints and simiilarity matrices for natural product databases.  

## Current status
- The analysis will be run on the COCONUT database (https://coconut.naturalproducts.net/). After removing compounds that cannot be parsed by RDKIT  
and the ones lacking taxonomy information, there are approx 120k NPs.  
- There are in total 5 NP taxonomy classes in this database:  
  - Plant (69%)    
  - Fungi (12%)   
  - Bacteria (10%)    
  - Marine (6%)  
  - Animal (0.8%)  
- There are currently 14 fingerprint types available:  
  - ECFP (rdkit)  
  - FCFP (rdkit)  
  - Avalon (rdkit)  
  - Topological torsion (rdkit)  
  - Atom-pair (rdkit)  
  - MACCS (rdkit)  
  - RDKIT (rdkit)  
  - Pubchem (cdk)  
  - CDK (cdk)  
  - ESTATE (cdk)  
  - Klekota-Roth (cdk)  
  - Lingo (cdk)  
  - MHFP (minhash)  
  - MAP4 (minhash)  
- The analysis will be done in batches of approx. 5000 compounds to deal with memory issues. MaxDiv sampling unfortunately requires choosing a FP type
to use in order to sample minimizing similarity (https://www.rdkit.org/docs/GettingStartedInPython.html#picking-diverse-molecules-using-fingerprints).  
As such, it would artifically bias the analysis against that fingerprint (if ECFP is used for MaxDiv, then it will have a bias towards lower average
tanimoto compared to the others, since the compounds were specifically picked to be diverse according to ECFP).  
To deal with this, I will sample randomly, but increase the number of sampling repetitions.  

## TODO
- [x] Add RDKIT fingerprints.  
- [x] Add CDK fingerprints.  
- [x] Wrap CDK fingerprints in Python.
- [x] Add minhash fingerprints.
- [x] Choose NP libraries to screen.  
- [x] Choose fingerprint comparison method.  
- [ ] Write fingerprint generation script.  
- [ ] Write similarity matrix generation script.  



