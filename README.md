# Estimating adaptive evolution across bacterial species 
This repository was created to upload the scripts used for the population genetics analysis for the paper:
"Recombination facilitates adaptive evolution in rhizobial soil bacteria", Cavassim et al., 2021.

## Primary steps (previously published in [Cavassim et al., 2020](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000351))
* Orthology: https://github.com/izabelcavassim/Rhizobium_analysis/#orthologous-identification-proteinortho
* Codon-aware alignment: https://github.com/izabelcavassim/Rhizobium_analysis/#codon-aware-alignment

## SNP-calling
python script: https://github.com/izabelcavassim/Rhizobium_analysis/#snp-calling

## Classifying synonymous and non-synonymous sites 

python script: hdf5_call_variants.py

Estimating transition transversion bias and non-synonysmous and synonymous sites
-----------------------
python script: ts_tv_calc.py


Computing the DFE and alpha
-----------------------
We used the method [GRAPES](https://github.com/BioPP/grapes) to estimate the DFE across pairs of species. 
We first computed the folded Site frequency spectrum of synonysmous and non-synonymous sites using the customized python script:

Grapes_files_pairwise_quantile.py

And then for each file input file generated we ran GRAPE as follow:



