#1 split the fasta files by genospecies
#2 run clonalframe (using the different species trees based on the concatenated genes)
#3 connect it to the rest of the pipeline
#4 code everything in the gwf manner
#5 update the github analysis
#6 write methodology in the draft
import pandas as pd
import numpy as np
import glob as glob
import os


def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip
	
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
		pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country, 'origin2':origin}
	return pop_map

def parse_fasta(filename, stratify=True, genospecies=None, pop_dict=None):
	
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	file_separe.remove('')
	parse_dict = {}
	header = []

	# Separating fasta files by genospecie
	if len(file_separe) == 196:

		for entry in file_separe:
			seq = entry.splitlines()
			header = seq[0] #these are the first elements of the list 

			strain_name = header.split("|")[1]
			#print strain_name
			gene_name = header.split("|")[0]

			# Modifying header with sequence position
			header = strain_name
			
			genome_id = seq[0].split('|')[1]

			seq = ''.join(seq[1:]) #joining the sequences 

			if stratify == True:
				if pop_dict[strain_name[0:4]]["genospecies"] == genospecies:
					parse_dict[header] = seq
			else:
				parse_dict[header] = seq
	return parse_dict

def write_fasta(tuple_fasta, directory):
	(fasta_name, dict_fasta) = tuple_fasta
	with open(directory+"{}".format(fasta_name), 'w') as f:
		for header, sequences in dict_fasta.items():
			f.write('>{}_{}\n'.format(fasta_name,header))
			chunks = [sequences[i:i+60] for i in range(0, len(sequences), 60)]
			for i in chunks:
				f.write('{}\n'.format(i))
		f.close()


pop_dict = parse_pop_map()

def split_fasta(genospecies, genospecies_dir, pop_dict=pop_dict):
	# Gene alignments: 
	genes_path = [(f) for f in glob.glob("/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/group_alns/*.fna")]
	for i in genes_path:
		gene_name = os.path.basename(i)
		gene_dict = parse_fasta(i, genospecies=genospecies, pop_dict=pop_dict)
		if len(gene_dict.keys()) > 0: 
			write_fasta((gene_name, gene_dict), directory=genospecies_dir)


split_fasta(pop_dict=pop_dict, genospecies="gsA", genospecies_dir="/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/")