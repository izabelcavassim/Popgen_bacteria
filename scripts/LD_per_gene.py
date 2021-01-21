# Linkage disequilibrium per gene and genospecies
# 1. Define the group of genes to be analysed
# 2. Segment it in genospecies
# 3. Calculate the LD for each genospecies
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
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

		#import numpy as np

file_candidates  = argv[1]  # #"/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/filtered_gene_statistics.csv"
metadata = argv[2] # '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'
results_dir = argv[3]  # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results/"
snp_file = argv[4] #'/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/newsnps_100.hdf5'

# syn_stats = pd.read_table(file_name_syn, delimiter=",")
# candidates = syn_stats["Gene"].tolist()


#Calculating LD just for chromosomal and chromids genes:
gene_groups = pd.read_csv(file_candidates)
candidates = gene_groups['Gene.group'].tolist()
print(len(candidates))


def parse_pop_map(file_name = metadata):
	#from itertools import izip   
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, origin, country, origin2 in zip(t['Seq ID'], t['Genospecies'], t['Country'], t['Origin2']):
		pop_map[str(strain_id)]={'genospecies':origin, 'country':country, 'origin':origin2}
	return pop_map


def calculating_avg(ld_dist_dict, distances, bin_size=100, gene=None, plot_figure=False, average_LD = True, candidates=candidates):

	# Initializing variables
	xs = []
	ys = []
	avg_r2s = list()
	plot_distances = list()
	avg_r2 = None

	#print(ld_dist_dict)
	for dist in distances:
		if ld_dist_dict[dist]['snp_count']>=1:
			avg_r2 = ld_dist_dict[dist]['r2_sum']/float(ld_dist_dict[dist]['snp_count'])
			avg_r2s.append(avg_r2)
			plot_distances.append(dist)

	plot_distances = sp.array(plot_distances)
	avg_r2s = sp.array(avg_r2s)
	print(avg_r2s)
	bins = sp.arange(0,max(plot_distances),bin_size)
	digitize = sp.digitize(plot_distances, bins) 

	print("These are the bins:")
	print(bins)
	print("These are the digitize")
	print(digitize)

	# Binning the results   
	for bin_i in range(len(bins)):
		bin_filter = digitize==(bin_i+1)
		xs.append(np.mean(plot_distances[bin_filter]))
		ys.append(np.nanmean(avg_r2s[bin_filter]))

	average_r2 = np.nanmean(ys)	
	print("This is the average r2")	
	print(average_r2)

	if average_LD == True:
		return(average_r2)


def LD_per_geno_per_gene(geno_species=[],
				 gt_hdf5_file= snp_file , candidates=candidates, min_maf=0.01, max_dist=1000):

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
			#print(gene)
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
	genes_accepted = []
	ld_decay_genes = []
	segregating_sites = []


	#for gene in gene_groups[0:100]:
	for i, gene in enumerate(gene_groups):
		if i % 100 == 0:
			print('Working on gene nr. %d' % i)
		char_name = "group" + gene
		gene = str(gene)

		if char_name in candidates:

			distances = range(1,max_dist)
			ld_dist_dict = {}
			for dist in distances:
				ld_dist_dict[dist]={'r2_sum':0.0, 'snp_count':0.0}


			print('Working on gene group: %s'%gene)
			# Gene strain list
			strains_list = ag[gene]['strains'][...]

			# Looking at specific genospecies
			gs_list = []
			for strain in strains_list:
				#strain = str(strain)
				strain=str(strain,'utf-8')
				gs_list.append((pop[strain]['genospecies'])) # +pop[strain]['origin'])

			# Transforming the strain list in array
			strains_list = np.asarray(strains_list)
			gs_filter1, gs_filter2 = [sp.in1d(gs_list,[gs]) for gs in geno_species]
			
			# Extracting the nucleotide sequences
			g1 = ag[gene]
			codon_snps = g1['codon_snps'][...]
			codon_snps = codon_snps[:,gs_filter1]
			norm_codon_snps = sp.transpose(codon_snps)
			#print(norm_codon_snps.shape)

			print("filtering low MAF snps")
			freqs = sp.mean(norm_codon_snps,0)
			norm_codon_snps = (norm_codon_snps-freqs)/sp.sqrt(freqs*(1-freqs))
			norm_codon_snps = sp.transpose(norm_codon_snps)
			mafs = sp.minimum(freqs,1-freqs)
			maf_filter = mafs>min_maf

			print(" Number of SNPs")
			print(norm_codon_snps.shape)

			print("Calculating LD")
			# Inposing a threshold of number of SNPs
			#if sp.sum(maf_filter)>=0:
			if norm_codon_snps.shape[0] >= 10: # at least 10 SNPs
				print(norm_codon_snps.shape)

				all_norm_snps = norm_codon_snps
				all_positions = g1['codon_snp_positions'][...]
				norm_snps = all_norm_snps[maf_filter]
				positions = all_positions[maf_filter]
				M,N = norm_snps.shape
				#print("LD shape")
				print(norm_snps.shape)
				is_synonimous_snp = g1['is_synonimous_snp'][...]
				is_nonsynonimous_snp = ~is_synonimous_snp
				syn_snp_filter = is_synonimous_snp*maf_filter
				nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter	
					
				if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
					all_norm_snps = norm_codon_snps
					all_positions = g1['codon_snp_positions'][...]
					norm_snps = all_norm_snps[maf_filter]
					positions = all_positions[maf_filter]
					M,N = norm_snps.shape
								
					ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
					#print(ld_mat.shape)
					assert M==len(positions), 'A bug detected.'
					for i in range(M-1):
						for j in range(i+1,M):
							dist = positions[j] - positions[i]
							if dist<max_dist:
								#print(ld_mat[i,j]**2)
								ld_dist_dict[dist]['r2_sum']+=ld_mat[i,j]**2
								ld_dist_dict[dist]['snp_count']+=1.0
		else:
			pass
		dct_sum = {k: sum(v.values()) for k, v in ld_dist_dict.items()}
		if sum(dct_sum.values()) > 0:
			genes_accepted.append(char_name)
			t = calculating_avg(ld_dist_dict, distances, gene=char_name)
			ld_decay_genes.append(t)
			segregating_sites.append(M)
			#else:
			#	continue
	dict_df = {"Genes":genes_accepted, "LD":ld_decay_genes, "Segregating_sites":segregating_sites}
	df = pd.DataFrame(dict_df)
	df.to_csv(results_dir+"LD_decay_{}_average_r2_10_snps_min_maf_0.01.csv".format(geno_species[0]))
	return(df)

# Estimating the LD per genospecies (gsA, gsB, gsC, gsD, and gsE)
# Analysing of all species
LD_per_geno_per_gene(["gsA", "gsA"])
LD_per_geno_per_gene(["gsB", "gsB"])
LD_per_geno_per_gene(["gsC", "gsC"])
LD_per_geno_per_gene(["gsD", "gsD"])
LD_per_geno_per_gene(["gsE", "gsE"])
