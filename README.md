# Popgen_bacteria
This repository was created to upload the scripts used for population genetics analysis of the Population genetics paper of sympatric species of Rhizobium.

## Primary steps 
* Orthology: https://github.com/izabelcavassim/Rhizobium_analysis/#orthologous-identification-proteinortho
* Codon-aware alignment: https://github.com/izabelcavassim/Rhizobium_analysis/#codon-aware-alignment

## SNP-calling
python script: https://github.com/izabelcavassim/Rhizobium_analysis/#snp-calling

## Classifying synonymous and non-synonymous sites 

python script: hdf5_call_variants.py

Estimating transition transversion bias
-----------------------
python script: ts_tv_calc.py


Estimating non-synonysmous and synonymous sites 
-----------------------

``` bash
for file in /home/mica16/NChain/faststorage/rhizobium/SFS_data/group_alns/{.,}*; do 
	basename=${file##*/}
	cp $file . 
	# replacing "-"'s (indels) by X
	python -c "import sys;lines=sys.stdin.read();print lines.replace('-','X')" < $basename > 'new_'$basename |

	# run the software to calculate the synonymous and non-synonymous sites
	qx --no-scratch -c 4 -m 8g -t 00:50:00 ./ratioSite 'new_'$basename 1
done
```
