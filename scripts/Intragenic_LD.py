#

# -*- coding: utf-8 -*-
# Site Frequency Spectrum
# 1. Define the group of genes to be analysed
# 2. Segment it in genospecies
# 3. Calculate the SFS for each genospecies
import h5py 
import pandas as pd
import numpy as np
import scipy as sp
from sys import argv
#import matplotlib
from collections import Counter
#from matplotlib import rcParams
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Tahoma']
#matplotlib.rcParams.update({'font.size': 10})
#import matplotlib.pyplot as plt
from decimal import *


filtered_genes = argv[1] #"/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/filtered_gene_statistics.csv"
metadata = argv[2] #'/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'
results_dir = argv[3] #'/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results'
snp_file = argv[4] #'/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/newsnps_100.hdf5'
maf = argv[5]
bin_size = argv[6]

def parse_pop_map(file_name = metadata):
	#from itertools import izip
	
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, sara_id, origin, country in zip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
		pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country, 'origin2':origin}
	return pop_map

def gen_ld_plots(snps_hdf5_file = snp_file , 
				 max_dist=2000, min_maf=0, bin_size=10,
				 fig_dir = results_dir, filter_pop=None, genes=filtered_genes):

	#Calculating LD just for chromosomal and chromids genes:
	gene_groups = pd.read_csv(genes)
	chrom_genes = gene_groups['Gene.group'].tolist()

	pop_map = parse_pop_map()

	xs = []
	ys = []

	#from itertools import izip
	h5f = h5py.File(snps_hdf5_file, mode="r")
	gene_groups = sorted(h5f.keys())
	ld_dist_dict = {'all':{}, 'nonsyn':{}, 'syn':{}}
	distances = range(0,max_dist)
	
	for dist in distances:
		ld_dist_dict['all'][dist]={'r2_sum':0.0, 'snp_count':0.0}
		#ld_dist_dict['nonsyn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
		#ld_dist_dict['syn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
	
	for i, gg in enumerate(chrom_genes):
	#for i, gg in enumerate(gene_groups):
		gg = gg[5::]
		print(gg)
		if gg in chrom_genes:
			print(gg)
		#gg = str(gg.encode('utf-8'))
		#print(type(gg))
		if i%100==0:
			print('%d: Gene %s'%(i,gg))

		g = h5f[str(gg)] 
		print(g)
				
		# Look at genes that have at least 10 SNPS
		if g['codon_snp_freqs'].size>10:

			if filter_pop is not None:
				strains = g['strains'][...]

				indiv_filter = sp.zeros((len(strains)),dtype='bool8')

				for s_i, s in enumerate(strains):
					try:
						s = str(s,'utf-8')
						if pop_map[s]['genospecies']==filter_pop:								
							indiv_filter[s_i]=True
					except:
						continue
				if sp.sum(indiv_filter)<2:
					continue
			
					
				codon_snps = g['codon_snps'][...]
				print(codon_snps)

				codon_snps = codon_snps[:,indiv_filter]
				print(codon_snps.shape)
				norm_codon_snps = sp.transpose(codon_snps)
				freqs = sp.mean(norm_codon_snps,0)
				norm_codon_snps = (norm_codon_snps-freqs)/sp.sqrt(freqs*(1-freqs))
				norm_codon_snps = sp.transpose(norm_codon_snps)
				mafs = sp.minimum(freqs,1-freqs)
				maf_filter = mafs>min_maf
				if sp.sum(maf_filter)>1:

					all_norm_snps = norm_codon_snps
					all_positions = g['codon_snp_positions'][...]
					norm_snps = all_norm_snps[maf_filter]
					positions = all_positions[maf_filter]
					M,N = norm_snps.shape
					is_synonimous_snp = g['is_synonimous_snp'][...]
					is_nonsynonimous_snp = ~is_synonimous_snp
					syn_snp_filter = is_synonimous_snp*maf_filter
					nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter
			
					if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
						all_norm_snps = norm_codon_snps
						all_positions = g['codon_snp_positions'][...]
						norm_snps = all_norm_snps[maf_filter]
						positions = all_positions[maf_filter]
						M,N = norm_snps.shape
								
						ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
						assert M==len(positions), 'A bug detected.'
						for i in range(M-1):
							for j in range(i+1,M):
								dist = positions[j] - positions[i]
								if dist<max_dist:
									ld_dist_dict['all'][dist]['r2_sum']+=ld_mat[i,j]**2
									ld_dist_dict['all'][dist]['snp_count']+=1.0

	print(ld_dist_dict)
	pairs = 0
	#for plot_type in ld_dist_dict.keys():
	avg_r2s = []
	plot_distances = []
	for dist in distances:
		if ld_dist_dict['all'][dist]['snp_count']>=1:
			avg_r2 = ld_dist_dict['all'][dist]['r2_sum']/float(ld_dist_dict['all'][dist]['snp_count'])
			pairs += 1
			avg_r2s.append(avg_r2)
			plot_distances.append(dist)
			
	plot_distances = sp.array(plot_distances)
	avg_r2s = sp.array(avg_r2s)

	print(avg_r2s)
	bins = sp.arange(0,max(plot_distances),bin_size)
	digitize = sp.digitize(plot_distances, bins)    
	for bin_i in range(len(bins)):
		bin_filter = digitize==(bin_i+1)
		if len(plot_distances[bin_filter])>0:
			xs.append(sp.mean(plot_distances[bin_filter]))
			ys.append(sp.mean(avg_r2s[bin_filter]))

		# plt.plot(xs, ys, color='k', linestyle='None', marker='.', alpha=0.5)
		# plt.xlabel(r'Pairwise distance ($d$)')
		# plt.ylabel(r'Squared correlation ($r^2$)')
		# if filter_pop is not None:
		# 	plt.title('LD decay of 0.99 < ANI <= 01')
		# 	plt.savefig('%s/ld_%s_codons_nuc_0.99_1_gsA_chromosome_maf_01_core_%s.pdf'%(fig_dir,plot_type,filter_pop))

	plot_list = pd.DataFrame(
	{'X': xs,
	 'Y': ys,
	})
	plot_list.to_csv("{dir_res}/plotting_intergenic_LD_{maf}_{bin_size}_{geno}.csv".format(dir_res=fig_dir, maf=min_maf, bin_size=bin_size, geno=filter_pop))
	return(plot_list)

gsA = gen_ld_plots(filter_pop='gsA')
gsB = gen_ld_plots(filter_pop='gsB')
gsC = gen_ld_plots(filter_pop='gsC')
gsD = gen_ld_plots(filter_pop='gsD')
gsE = gen_ld_plots(filter_pop='gsE')
