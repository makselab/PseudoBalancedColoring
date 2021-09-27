# SymmetryGroupFactorization

.py files contain classes that work with symmetry groups using NAUTY, Sage and SymPy.

"MIPoutputAnalyzer.ipynb" notebook contains some examples of how this analysis can be performed.
If you are working with the output from the MIP use the second cell. If you are working with the output from the MIP file by file use the third cell. If you are working with the edge file, use the fourth cell.
This notebook needs to be ran from Sage like so "sage -n jupyter". Make sure that PyNauty and SymPy are installed and included in the Sage environment.

"pythonIndicesAnalysis.R" can be used to reproduce plots in the paper.
Run it on the folder with the outputs from MIP. It will collect the data, add Kamei colors to the node files and make the plot.

"manualRepairIndices" contains the indices obtained on the gap junction circuits that Morone and Makse used in their Nature Communications paper.

"test" folder contains some sample files that can be worked with.
