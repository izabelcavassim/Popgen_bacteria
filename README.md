# Popgen_bacteria
This repository was created to upload the scripts used for population genetics analysis of the XX paper

## SNP-calling
python script: 

## Classifying synonymous and non-synonymous sites 

python script: hdf5_call_variants.py


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
