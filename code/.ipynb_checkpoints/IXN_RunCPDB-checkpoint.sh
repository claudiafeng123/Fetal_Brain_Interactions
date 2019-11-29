cd ~/rotation1/

## ---- SetVariables

dataFolder="data/IXN/CPDB/cpdb_input/"
outFolder="data/IXN/CPDB/cpdb_results/"
figureFolder="figures/IXN/"
toolsFolder="../tools/"

## ---- SetUp

#python -m venv $toolsFolder
#source ${toolsFolder}cpdb-venv/bin/activate
#pip install cellphonedb

#pip install cffi==1.12.3

## ---- RunCPDB

#activate CPDB virtual environment
source ${toolsFolder}cpdb-venv/bin/activate

#for testing, use:\
#i=${dataFolder}fetus_meta.txt
#i=${dataFolder}Giandomenico2019_H9-75d-Batch1_meta.txt
for i in $dataFolder*_meta.txt
do
    ID=`echo $i | awk -F "/" '{print $5}' | sed 's/_meta.txt//'`
    echo $ID

    #run cellphonedb
    cellphonedb method statistical_analysis $dataFolder${ID}_meta.txt $dataFolder${ID}_counts.txt --output-path ${outFolder} --project-name $ID --threads=16 --pvalue=0.05 --threshold=0.1
    
    #dotplot
    cellphonedb plot dot_plot --output-path ${figureFolder} --means-path ${outFolder}$ID/means.txt --pvalues-path ${outFolder}$ID/pvalues.txt --output-name CPDB_dot-plot_${ID}.pdf

    #heatmap
    cellphonedb plot heatmap_plot $dataFolder${ID}_meta.txt --output-path ${figureFolder} --pvalues-path ${outFolder}$ID/pvalues.txt --count-name CPDB_heatmap_${ID}.pdf

done


