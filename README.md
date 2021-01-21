# Estimating adaptive evolution among bacterial species 
This repository was created to upload the scripts used for the population genetics analyses of the paper:
["Recombination facilitates adaptive evolution in rhizobial soil bacteria", Cavassim et al., 2021.](https://doi.org/10.1101/2021.01.20.427438)


## Primary steps (previously published in [Cavassim et al., 2020](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000351))
* Orthology search: https://github.com/izabelcavassim/Rhizobium_analysis/#orthologous-identification-proteinortho
* Codon-aware alignment of orthologous genes: https://github.com/izabelcavassim/Rhizobium_analysis/#codon-aware-alignment

## SNP-calling
python script: https://github.com/izabelcavassim/Rhizobium_analysis/#snp-calling

## Classifying synonymous and non-synonymous sites 

python script: [hdf5_call_variants.py](https://github.com/izabelcavassim/Popgen_bacteria/blob/master/scripts/hdf5_call_variants.py)

Estimating transition transversion bias and non-synonysmous and synonymous sites
-----------------------
python script: [ts_tv_calc.py](https://github.com/izabelcavassim/Popgen_bacteria/blob/master/scripts/ts_tv_calc.py)


Computing the DFE and alpha
-----------------------
We used the method [GRAPES](https://github.com/BioPP/grapes) to estimate the DFE across pairs of species. 
We first computed the folded Site frequency spectrum of synonysmous and non-synonymous sites using the customized python script:

Grapes_files_pairwise_quantile.py

To produce a text file ({species1}_{species2}_grapes.txt) that looks like this one:

``` R
gsE+gsD (548.000000 genes)
'all_genes'	11	547966.6905329467	365	257	354	107	259	149830.30946705345	2698	2220	2469	1062	1952	547966.6905329467	963	149830.30946705345	8810
```
In which species_1 (this case gsE) is the focal group (polymorphism) and species 2 (in this case gsD) is the outgroup (divergence). 
A description of each entry is found [here](https://github.com/BioPP/grapes#example-input-files-for-grapes). 

And then for each input file generated we ran [GRAPES](https://github.com/BioPP/grapes) as follow:

``` bash
/home/mica16/grapes-izabel-work/bin/grapes -in {inputs_dir}{species1}_{species2}_grapes.txt -out {results_dir}{species1}_{species2}_grapes_output.txt -nb_rand_start 20 -model all
```
The same procedure is done for computing alpha across recombination classes. 

Workflow
-----------------------
For my own sake I created a python workflow that combines all the described analyses above. I used the software [gwf](https://gwf.app/) to build and run this workflow.

Workflow pipeline (workflow.py), scripts and files connected to it are found [here](https://github.com/izabelcavassim/Popgen_bacteria/blob/master/scripts/)

Feel free to contact me if you want to apply these scripts in your own dataset or if you have any questions. 
