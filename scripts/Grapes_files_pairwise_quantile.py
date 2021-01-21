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
from sys import argv
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
matplotlib.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
from decimal import *


filtered_genes = argv[1] # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/filtered_gene_statistics.csv"
filtered_name_syn  = argv[2]  # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results/Synonymous_sites_counts.csv"
filtered_name_non_syn  = argv[3]
metadata = argv[4] # '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'
results_dir = argv[5]  # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results/"
snp_file = argv[6] #'/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/newsnps_100.hdf5'
quantiles = argv[7]
input_recombination_files = argv[8]
print(quantiles)
def parse_pop_map(file_name = metadata):
	#from itertools import izip   
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
	from collections import OrderedDict

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
	SFS_counts = dict(sorted(SFS_counts.items()))
	return SFS_counts

#Filtered genes
gene_groups = pd.read_csv(filtered_genes)
candidates = gene_groups['Gene.group'].tolist()

def divergence(geno_species=[], bin_size=0.2,
				 gt_hdf5_file= snp_file, candidates=candidates, filename="test.txt", syn_file = filtered_name_syn, non_syn_file = filtered_name_non_syn, results_dir=results_dir):


	file_name_syn = filtered_name_syn
	syn_stats = pd.read_table(file_name_syn, delimiter=",")
	syn_stats = syn_stats.set_index('Gene')["Synonymous_sites"].to_dict()

	# Non-Synonymous sites counts
	file_name_nonsyn = filtered_name_non_syn
	nonsyn_stats = pd.read_table(file_name_nonsyn, delimiter=",")
	nonsyn_stats = nonsyn_stats.set_index('Gene')["Non_Synonymous_sites"].to_dict()


	non_syn_sites = 0
	syn_sites = 0

	# The hdf5 contains:
	#[u'aacids', u'blosum62_scores', u'codon_snp_freqs', u'codon_snp_positions', u'codon_snps', u'codons', u'diversity',
	# u'dn_ds_ratio', u'freqs', u'is_synonimous_snp', u'norm_codon_snps', u'norm_snps', u'nts', u'num_non_syn_sites', 
	# u'num_syn_sites', u'num_vars', u'raw_snp_positions', u'raw_snps', u'snp_positions', u'snps', u'strains']

	pop = parse_pop_map()
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file, mode="r")

	ag = h5f
	#print(ag.keys())

	gene_big_groups = sorted(ag.keys())
	gene_groups = list()

	# Taking just the core genes
	for gene in gene_big_groups:
		if len(ag[gene]['strains']) == 196:
			gene_groups.append(gene)
	#print(gene_groups)
	
	# Names
	# sorted strains by genospecies
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])

	# Deleting some strains
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')
	strains_list = strains_names


	# Preparing variables 
	g1_list = []
	g1_syn = []
	g1_non = []
	g2_list = []
	g2_syn = []
	g2_non = []
	divergence_syn = 0
	divergence_nonsyn = 0
	divergence_syn_sites = 0
	divergence_nonsyn_sites = 0
	shared_pol_syn = 0
	shared_pol_nonsyn = 0
	combined = list()
	counts = 0

	decoder = np.vectorize(lambda x: x.decode("utf-8"))
	print(gene)
	for i, gene in enumerate(gene_groups):
		if i % 100 == 0:
			print('Working on gene nr. %d' % i)
		char_name = "group" + gene

		## Restricting our data to only specific candidates
		if char_name in candidates:
			counts += 1
			strains_list = ag[gene]['strains'][...]
			strains_list = decoder(strains_list)

			# MODIFYYYYYYYYYYY THIS PART (I wanna build a synonymous and non-synonymous files with all gene counts )
			# Adding synonymous non-synonymous sites (possible sites)
			try:
				char_name_mod = char_name
				#print(char_name_mod)
				#print(nonsyn_stats[char_name_mod])
				non_syn_sites += nonsyn_stats[char_name_mod]
				syn_sites += syn_stats[char_name_mod]
				#print(syn_sites)
			except:
				continue

			# Looking at specific genospecies
			gs_list = []
			for strain in strains_list:
				gs_list.append((pop[strain]['genospecies'])) 

			# Transforming the strain list in array
			strains_list = np.asarray(strains_list)
			gs_filter1, gs_filter2 = [sp.in1d(gs_list,[gs]) for gs in geno_species]
		
			############## Extracting species indexes ################### 		
			gs1 = strains_list[gs_filter1]
			gs2 = strains_list[gs_filter2]
			total_gs = np.append(gs1,gs2)

			############## Extracting the nucleotide sequences ################### 
			print(gene)
			g1 = ag[gene]
			bla = g1['nts'][...].tolist()
			bla = list(itertools.chain.from_iterable(bla))
			#print(b" ".join(bla))
			g1 = g1['codon_snps'][...].T
			g1bla = list(itertools.chain.from_iterable(g1))
			#print(g1bla)

			syn_index = ag[gene]['is_synonimous_snp'][...]
			print(syn_index)
			
			#### First species
			g1_geno = g1[gs_filter1, :]
			g1_list.append(g1_geno)

			#### Rows are strains and columns are SNPs ####
			g1_vector_syn = sum(g1_geno[:,syn_index]) # sum of minor (0) and major allele (1)
			g1_vector_nonsyn = sum(g1_geno[:,~syn_index])
			g1_syn.append(g1_geno[:,syn_index])
			g1_non.append(g1_geno[:,~syn_index])

			#### Second species
			g2_geno = g1[gs_filter2, :]
			g2_list.append(g2_geno)

			g2_vector_syn = sum(g2_geno[:,syn_index]) # sum of minor (0) and major allele (1)
			g2_vector_nonsyn = sum(g2_geno[:,~syn_index])
			g2_syn.append(g2_geno[:,syn_index])
			g2_non.append(g2_geno[:,~syn_index])


			############## Shared Polymorphisms ###################
			g1_max = g1_geno[:,syn_index].shape[0]
			g2_max = g2_geno[:,syn_index].shape[0]

			#### Synonymous
			# Find the places in the gene were the number of 1's or zeros is above 0 and below maximum (polymorphic sites):
			pol_sites_gs2_syn = np.where((g2_vector_syn < g2_max) & (g2_vector_syn > 0))
			pol_sites_gs1_syn = np.where((g1_vector_syn < g1_max) & (g1_vector_syn > 0))
			shared_pol_syn += len(np.intersect1d(pol_sites_gs1_syn, pol_sites_gs2_syn))

			#### Non-synonymous
			# Find the places in the gene were the number of 1's or zeros is above 0 and below maximum:
			pol_sites_gs2_nonsyn = np.where((g2_vector_nonsyn < g2_max) & (g2_vector_nonsyn > 0))
			pol_sites_gs1_nonsyn = np.where((g1_vector_nonsyn < g1_max) & (g1_vector_nonsyn > 0))
			shared_pol_nonsyn += len(np.intersect1d(pol_sites_gs1_nonsyn, pol_sites_gs2_nonsyn))

			############## Fixed differences ###################
			# Looking at divergence of species 1, 
			# Species 1 is complete divergent if species 1 is at its maximum and the other species is 0 
			# Or if is at its minimum and the other species is at its maximum

			#### Synonymous
			min_g1_syn = np.where(g1_vector_syn == 0)
			max_g1_syn = np.where(g1_vector_syn == g1_max)

			min_g2_syn = np.where(g2_vector_syn == 0)
			max_g2_syn = np.where(g2_vector_syn == g2_max)

			# Is at its maximum and the other species is at its minimum:
			divergence_syn += len(np.intersect1d(max_g1_syn,min_g2_syn))

			# Its minimum and the other species is at its maximum		
			divergence_syn += len(np.intersect1d(min_g1_syn,max_g2_syn))

			# Divergence sites: any place where both species are not zero or fixed
			zeros = len(np.intersect1d(min_g1_syn, min_g2_syn))
			ones = len(np.intersect1d(max_g1_syn, max_g2_syn))
			calc_syn = len(g1_vector_syn) - (zeros + ones)

			# Divergence sites: fixed in both species
			divergence_syn_sites += calc_syn 

			##### Non-synonymous
			min_g1_nonsyn = np.where(g1_vector_nonsyn == 0)
			max_g1_nonsyn = np.where(g1_vector_nonsyn == g1_max)

			min_g2_nonsyn = np.where(g2_vector_nonsyn == 0)
			max_g2_nonsyn = np.where(g2_vector_nonsyn == g2_max)

			# Is at its maximum and the other species is at its minimum:
			divergence_nonsyn += len(np.intersect1d(max_g1_nonsyn,min_g2_nonsyn))

			# Its minimum and the other species is at its maximum
			divergence_nonsyn += len(np.intersect1d(min_g1_nonsyn,max_g2_nonsyn))


	g1_conc = np.concatenate(g1_list, axis = 1) 
	g1_syn_conc = np.concatenate(g1_syn, axis = 1)
	g1_non_conc = np.concatenate(g1_non, axis = 1)

	print(g1_syn_conc.shape)
	all_sites = g1_syn_conc.shape[1] + g1_non_conc.shape[1]
	print(all_sites)

	print("Synonymous and non-synonymous sites")
	print(syn_sites, non_syn_sites)

	print('Number of genes analysed: %d' % counts)

	print(type(all_sites))
	print(geno_species)
	print('all sites (SNPs) %f' % all_sites)
	print('synonymous polymorphic sites %f' % syn_sites)
	print('non-synonymous polymorphic sites %f' % non_syn_sites)
	print('non-synonymous divergence sites %f' % non_syn_sites)
	print('non-synonymous divergence substitutions %f' % divergence_nonsyn)
	print('synonymous divergence sites %f' % syn_sites)
	print('synonymous divergence substitutions %f' % divergence_syn)
	print('shared synonymous polymoprhisms %f' % shared_pol_syn)
	print('shared nonsynonymous polymoprhisms %f' % shared_pol_nonsyn)

	sfs_codon =  num_segregating_sites(g1_conc)
	sfs_codon = pd.DataFrame.from_dict(sfs_codon, orient='index')

	print('SFS synonymous sites')
	sfs_syn =  num_segregating_sites(g1_syn_conc)
	sfs_syn_df = pd.DataFrame.from_dict(sfs_syn, orient='index')

	print('SFS non-synonymous sites')
	sfs_non =  num_segregating_sites(g1_non_conc)
	sfs_non_df = pd.DataFrame.from_dict(sfs_non, orient='index')

	# results combined
	combined.append("all_genes")
	combined.append(g1_max)
	combined.append(non_syn_sites)
	print(sfs_non)
	combined = combined + list(sfs_non.values())
	combined.append(syn_sites)
	print(sfs_syn)
	combined = combined + list(sfs_syn.values())
	combined.append(non_syn_sites)
	combined.append(divergence_nonsyn)
	combined.append(syn_sites)
	combined.append(divergence_syn)

	print(combined)

	# Creating the output for grapes:
	with open(results_dir+filename, 'w') as f:
		f.write("%s+%s (%f genes)\n" % (geno_species[0], geno_species[1], len(g1_list)))
		f.write("	".join(repr(e) for e in combined))
	f.close()
	
	df = pd.concat([sfs_syn_df, sfs_non_df], axis=1)
	df.columns = ['Syn', 'Non-syn']

	return([len(gene_groups), sfs_non, sfs_syn, shared_pol_syn, shared_pol_nonsyn])

def pairwise_files():
	# Analysing gsA
	genos = ['gsA', 'gsB','gsC', 'gsD', 'gsE']
	comb = list(itertools.permutations(genos, 2))

	results_shared_pol = list()
	for i in comb:
		lists_sfs = divergence(i, filename=i[0]+"_"+i[1]+'_'+"grapes.txt")
		results_shared_pol.append([i, lists_sfs[3], lists_sfs[4]])
	print(results_shared_pol)

def grapes_files_per_recombination(results=results_dir):

	recombination_quantiles = [(f) for f in glob.glob(input_recombination_files+"recombination_quantiles*")]


	genos = ['gsA', 'gsB','gsC', 'gsD', 'gsE']
	comb = list(itertools.permutations(genos, 2))

	for i in recombination_quantiles:

		#Recombination quantiles
		file_quantiles = i
		geno = i.split("/")[8][-7::]
		geno = geno[0:3]
		print(geno)

		rec_quantiles = pd.read_table(file_quantiles, delimiter=",")
			
		# First quantile:
		first = rec_quantiles.loc[rec_quantiles['quartile'] == 1]
		first['gene'] = first['gene'].astype(str)

		# Second quantile
		second = rec_quantiles.loc[rec_quantiles['quartile'] == 2]
		second['gene'] = second['gene'].astype(str)

		# Thrid quantile
		#if geno != "gsD":
		third = rec_quantiles.loc[rec_quantiles['quartile'] == 3]
		third['gene'] = third['gene'].astype(str)
		#fourth = rec_quantiles.loc[rec_quantiles['quartile'] == 4]
		#fourth['gene'] = fourth['gene'].astype(str)
		#print(first)

		# Only A and E for now
		for c in comb:
			if geno in c[0]:
			#print(c)
			#comb = list(itertools.permutations(genos, 2))
				divergence(c, filename=c[0]+"_"+c[1]+'_'+"_first_quantile_grapes.txt", candidates= first["gene"].tolist(), results_dir=results)
				divergence(c, filename=c[0]+"_"+c[1]+'_'+"_second_quantile_grapes.txt", candidates= second["gene"].tolist(), results_dir=results)
				divergence(c, filename=c[0]+"_"+c[1]+'_'+"_third_quantile_grapes.txt", candidates= third["gene"].tolist(), results_dir=results)
							#divergence(c, filename=c[0]+"_"+c[1]+'_'+"_fourth_quantile_grapes.txt", candidates= third["gene"].tolist(), results_dir=results)
				#if c[0] == "gsD" and geno in c:
							#divergence(c, filename=c[0]+"_"+c[1]+'_'+"_first_quantile_grapes.txt", candidates= first["gene"].tolist(), results_dir=results)
							#divergence(c, filename=c[0]+"_"+c[1]+'_'+"_second_quantile_grapes.txt", candidates= second["gene"].tolist(), results_dir=results)

if quantiles == "True":
	grapes_files_per_recombination()
elif quantiles == "False":
	pairwise_files()
