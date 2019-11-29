cd ~/rotation1/

## ---- SetVariables

cpdb_results="data/IXN/cpdb_results/"
figureFolder="figures/RSC/"
toolsFolder="../tools/"

LRPairs="data/RSC/LRPairs.txt"
IPColumns="data/RSC/IP_columns.txt"

## ---- DotPlots

#activate CPDB virtual environment
source ${toolsFolder}cpdb-venv/bin/activate

#for testing, use:\
ID=fetus
cellphonedb plot dot_plot --output-path ${figureFolder} --means-path ${cpdb_results}$ID/means.txt --pvalues-path ${cpdb_results}$ID/pvalues.txt --output-name CPDB_dot-plot_${ID}.pdf --rows ${LRPairs} --columns ${IPColumns}


