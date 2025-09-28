# GWAS_analysis

## Requirements
- NextFlow
- Python libraries
    - pandas
    - numpy
    - cyvcf2
    - matplotlib
    - seaborn
    - scikit-learn
    - scipy
    - scikit-allel
    - natsort
 
## Running
RUN in the folder with test data
<pre>
nextflow run main.nf
</pre>
    
## Results
Files are saved to a folder ./results
Description | Command
--- | ---
PCA.png | a scatter plot after PCA for genotypes and clustering by phenotypes 
tree.png | a phylogentic tree
gwas_pvalues.csv | P-values for GWAS analysis
manhattan_plot.png | Manhattan plot
QQ_plot.png | QQ plot
