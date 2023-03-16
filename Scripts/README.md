# Folder structure

- [FP_calc:](clean_coconut.csv) Contains all functions and supporting files to compute
all fingerprints. Currently the package can calculate:  
   - MACCS  
   - Extended connectivity  
   - Functional class  
   - RDKIT  
   - Atom Pair  
   - Topological torsion  
   - Avalon  
   - MinHashed  
   - MAP4 
   - PUBCHEM 
   - Daylight 
   - ESTATE 
   - Klekotha - Roth 
   - LINGO 
   - Depth First Search 
   - All Shortest Paths 
   - LSTAR 
   - Molprint 
   - Pharmacophoric pairs 
   - Pharmacophoric triplets  

- [cleanup.py:](cleanup.py) Contains utility functions used in `cleanup_script`.  

- [cleanup_script.py:](cleanup_script.py) Preprocesses the raw COCONUT database. For further
information, open the file or use the `--help` flag from the command line.  

- [clf.py:](clf.py) Contains utility functions used in `clf_script`.  

- [clf_script.py:](clf_script.py) Calculates classification metrics for all classes
using all precomputed fingerprints. For further information, open the file or use the `--help` flag from the command line.  

- [fp_script.py:](fp_script.py) Computes the chosen fingerprint set. For further information, open the file or use the `--help` flag from the command line.  

- [sim_search.py:](sim_search.py) Contains utility functions used in `sim_search_script`.  

- [sim_search_script.py:](sim_search_script.py) Runs the similarity search using all precomputed fingerprints. For further information,
 open the file or use the `--help` flag from the command line.  

- [utils.py:](utils.py) Contains utility functions used throughout the repository.  

