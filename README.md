# SymmetryGroupFactorization

Authors:
- Decomposition code: **Ian Leifer, Levich Institute and Physics Department, City College of New York, New York, NY 10031**
- IP code: **David Phillips, Mathematics Department, United States Naval Academy, Annapolis, MD, 21402**

Advisor: **Hernan Makse, Levich Institute and Physics Department, City College of New York, New York, NY 10031**

.py files contain classes that work with symmetry groups using NAUTY, Sage and SymPy.

"MIPoutputAnalyzer.ipynb" notebook contains some examples of how this analysis can be performed.
If you are working with the output from the MIP use the second cell. If you are working with the output from the MIP file by file use the third cell. If you are working with the edge file, use the fourth cell.
This notebook needs to be ran from Sage like so "sage -n jupyter". Make sure that PyNauty and SymPy are installed and included in the Sage environment.

"pythonIndicesAnalysis.R" can be used to reproduce plots in the paper.
Run it on the folder with the outputs from MIP. It will collect the data, add Kamei colors to the node files and make the plot.

"manualRepairIndices" contains the indices obtained on the gap junction circuits that Morone and Makse used in their Nature Communications paper.

"test" folder contains some sample files that can be worked with.
