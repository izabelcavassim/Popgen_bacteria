# Site Frequency Spectrum
# 1. Define the group of genes to be analysed
# 2. Segment it in genospecies
# 3. Calculate the SFS for each genospecies
import h5py 
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib
import itertools
import glob
from fractions import *
from collections import Counter
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
matplotlib.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
from decimal import *
from sys import argv

file_candidates  = argv[1]  # #"/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/filtered_gene_statistics.csv"
metadata = argv[2] # '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'
results_dir = argv[3]  # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results/"
snp_file = argv[4] #'/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/newsnps_100.hdf5'

#Calculating LD just for chromosomal and chromids genes:
gene_groups = pd.read_csv(file_candidates)
candidates = gene_groups['Gene.group'].tolist()


def parse_pop_map(file_name = metadata):
	#from itertools import zip   
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, origin, country, origin2 in zip(t['Seq ID'], t['Genospecies'], t['Country'], t['Origin2']):
		pop_map[str(strain_id)]={'genospecies':origin, 'country':country, 'origin':origin2}
	return pop_map

def num_segregating_sites(gene_matrix):
	"""
	Input snp matrix
	Returns the raw number of segregating sites (polymorphic sites).
	Sum over collumns, if sum != 0 or sum != nrow(matrix) : segregating site
	"""

	#for i in len(gene_matrix.shape[0] - 1):
	#gene_matrix = numpy.delete(gene_matrix, (0), axis=0)

	freqs = sp.mean(gene_matrix, 0)
	mafs = sp.minimum(freqs, 1 - freqs)

	# Filtering very rare variants?
	# maf_filter = mafs > 0.001
	# mafs = mafs[maf_filter]

	sum_list = mafs * gene_matrix.shape[0]
	data = [float(Decimal("%.2f" % e)) for e in sum_list]

	SFS = Counter(data)

	del SFS[0.0] # all fixed values

	total = sum(SFS.values(), 0.0)
	
	SFS_freq = {k: v / total for k, v in SFS.items()}
	SFS_counts = {k: v for k, v in SFS.items()}

	print(SFS_counts)
	return SFS_counts


def SFS_per_geno(geno_species=[], bin_size=0.2,
				 gt_hdf5_file= snp_file , candidates=candidates):

	# The hdf5 contains:

	#[u'aacids', u'blosum62_scores', u'codon_snp_freqs', u'codon_snp_positions', u'codon_snps', u'codons', u'diversity',
	# u'dn_ds_ratio', u'freqs', u'is_synonimous_snp', u'norm_codon_snps', u'norm_snps', u'nts', u'num_non_syn_sites', 
	# u'num_syn_sites', u'num_vars', u'raw_snp_positions', u'raw_snps', u'snp_positions', u'snps', u'strains']

	pop = parse_pop_map()
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file, mode="r")

	ag = h5f

	gene_big_groups = sorted(ag.keys())
	gene_groups = list()

	# Taking just the core genes
	for gene in gene_big_groups:
		#print len(ag[gene]['strains'])
		if len(ag[gene]['strains']) == 196:
			gene_groups.append(gene)

	print('Number of genes analysed: %f' % len(gene_groups))
	# Names
	# sorted strains by genospecies

	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	#print 'These are the strains evaluated', strains_names
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')
	strains_list = strains_names

	g1_list = []
	g1_syn = []
	g1_non = []

	count = 0
	#for gene in gene_groups[0:100]:
	for i, gene in enumerate(gene_groups):
		if i % 100 == 0:
			print('Working on gene nr. %d' % i)
		char_name = "group" + gene
		if char_name in candidates:
			count += 1

			#print 'Working on gene group: %s'%gene
			# Gene strain list
			strains_list = ag[gene]['strains'][...]

			# Looking at specific genospecies
			gs_list = []
			for strain in strains_list:
				strain=str(strain,'utf-8')
				gs_list.append((pop[strain]['genospecies'])) # +pop[strain]['origin'])

			# Transforming the strain list in array
			strains_list = np.asarray(strains_list)
			gs_filter1, gs_filter2 = [sp.in1d(gs_list,[gs]) for gs in geno_species]
				
			gs1 = strains_list[gs_filter1]
			gs2 = strains_list[gs_filter2]
			total_gs = np.append(gs1,gs2)

			# Extracting the nucleotide sequences
			g1 = ag[gene]
			g1 = g1['codon_snps'][...].T

			try:
				g1_geno = g1[gs_filter1, :]
				g1_list.append(g1_geno)

				syn_index = ag[gene]['is_synonimous_snp'][...]
				g1_syn.append(g1_geno[:,syn_index])
				g1_non.append(g1_geno[:,~syn_index])
			except:
				pass

	print(len(g1_list))
	print(count)
	g1_conc = np.concatenate(g1_list, axis = 1) 
	g1_syn_conc = np.concatenate(g1_syn, axis = 1)
	g1_non_conc = np.concatenate(g1_non, axis = 1)
	
	print('all %f' % g1_conc.shape[1])	
	print('synonymous %f' % g1_syn_conc.shape[1])
	print('non-synonymous %f' % g1_non_conc.shape[1])

	sfs_codon =  num_segregating_sites(g1_conc)
	sfs_codon = pd.DataFrame.from_dict(sfs_codon, orient='index')

	sfs_syn =  num_segregating_sites(g1_syn_conc)
	sfs_syn = pd.DataFrame.from_dict(sfs_syn, orient='index')

	sfs_non =  num_segregating_sites(g1_non_conc)
	sfs_non = pd.DataFrame.from_dict(sfs_non, orient='index')
	
	df = pd.concat([sfs_syn, sfs_non], axis=1)
	df.columns = ['Syn', 'Non-syn']


	df.to_csv(results_dir+"SFS_{}_counts.csv".format(geno_species[0]))
	t = 'SFS for genospecies %s organic field' % geno_species[0]
	ax = df.plot(kind='bar',  title = t)
	ax.set_xlabel('Minor allele count')
	ax.set_ylabel('Fraction of SNPs')
	ax.set_title('SFS of %s' % geno_species[0])
	plt.tight_layout()
	plt.savefig(results_dir+'SFS_%s.pdf' % geno_species[0])
	#plt.show()
	return([len(gene_groups), sfs_non, sfs_syn])

# Estimating the folded Site Frequency Spectrum 
# Analysing gsA
SFS_per_geno(["gsA", "gsA"])

# Analysing gsB
SFS_per_geno(["gsB", "gsB"])

# Analysing gsC
SFS_per_geno(["gsC", "gsC"])

# Analysing gsD
SFS_per_geno(["gsD", "gsD"])

# Analysing gsE
SFS_per_geno(["gsE", "gsE"])
