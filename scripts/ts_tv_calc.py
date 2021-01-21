# Calculating transition transversion ratio
# Input alignment
# Calculate pairwise
# Take the average
# Try to calculate with species

from itertools import combinations
import glob
import os
import pandas as pd
from collections import OrderedDict 
from Bio.codonalign.codonalphabet import default_codon_table
from math import log
import numpy as np
import random
import matplotlib.pyplot as plt
from sys import argv

filtered_genes_directory  = argv[1]  # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/filtered_gene_statistics.csv"
group_alignment_dir = argv[2] #"/Users/PM/Dropbox/Cavassim_et_al_2019_Rhizobium_data/group_alns/"
results_dir = argv[3] = # "/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/Results/"

def parse_fasta(filename, stratify=False, genospecies=None, pop_dict=None, basename=None):
		
		file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
		file_separe = file.split('>') #spliting each entry by the > 
		file_separe.remove('')
		parse_dict = OrderedDict()
		header = []


		# Separating fasta files by genospecie
		if len(file_separe) <= 196:

				for entry in file_separe:
						seq = entry.splitlines()
						header = seq[0] #these are the first elements of the list 

						strain_name = header.split("|")[1]
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

		# Writing new version of the nucleotide sequence:
		# file = open("/Users/PM/Dropbox/PHD/Pop_gen_rhizobium_paper2/{gene}.fasta".format(gene=basename), "w")
		# for header, seq in parse_dict.items():
		# 	file.write(">" + str(header) + "\n" + str(seq) + "\n")
		return parse_dict

def _get_pi(seq1, seq2, cmethod, codon_table=default_codon_table): 
	"""Obtain codon frequency dict (pi) from two codon list (PRIVATE). 
  
	  This function is designed for ML method. Available counting methods 
	 (cfreq) are F1x4, F3x4 and F64. 
	""" 
	# TODO: 
	# Stop codon should not be allowed according to Yang. 
	# Try to modify this! 
	pi = {} 
	if cmethod == "F1x4": 
		fcodon = {"A": 0, "G": 0, "C": 0, "T": 0} 
		for i in seq1 + seq2: 
			if i != "---": 
				for c in i: 
					fcodon[c] += 1 
		tot = sum(fcodon.values()) 
		fcodon = {j: k / tot for j, k in fcodon.items()} 
		for i in codon_table.forward_table.keys() + codon_table.stop_codons: 
			if "U" not in i: 
				pi[i] = fcodon[i[0]] * fcodon[i[1]] * fcodon[i[2]] 
	elif cmethod == "F3x4": 
		# three codon position 
		fcodon = [{"A": 0, "G": 0, "C": 0, "T": 0}, {"A": 0, "G": 0, "C": 0, "T": 0}, {"A": 0, "G": 0, "C": 0, "T": 0}] 
		for i in seq1 + seq2: 
			if i != "---": 
				fcodon[0][i[0]] += 1 
				fcodon[1][i[1]] += 1 
				fcodon[2][i[2]] += 1 
		for i in range(3): 
			tot = sum(fcodon[i].values()) 
			fcodon[i] = {j: k / tot for j, k in fcodon[i].items()} 
		for i in list(codon_table.forward_table.keys()) +  codon_table.stop_codons: 
			if "U" not in i: 
				pi[i] = fcodon[0][i[0]] * fcodon[1][i[1]] * fcodon[2][i[2]] 
	elif cmethod == "F61": 
		for i in codon_table.forward_table.keys() + codon_table.stop_codons: 
			if "U" not in i: 
				pi[i] = 0.1 
		for i in seq1 + seq2: 
			if i != "---": 
				pi[i] += 1 
		tot = sum(pi.values()) 
		pi = {j: k / tot for j, k in pi.items()} 
	return pi 

def _count_site_YN00(codon_lst1, codon_lst2, pi, k, codon_table=default_codon_table):
	"""Site counting method from Ina / Yang and Nielsen (PRIVATE).
	Method from `Ina (1995)`_ as modified by `Yang and Nielsen (2000)`_.
	This will return the total number of synonymous and nonsynonymous sites
	and base frequencies in each category. The function is equivalent to
	the ``CountSites()`` function in ``yn00.c`` of PAML.
	.. _`Ina (1995)`: https://doi.org/10.1007/BF00167113
	.. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236
	"""
	if len(codon_lst1) != len(codon_lst2):
		raise RuntimeError(
			"Length of two codon_lst should be the same (%d and %d detected)"
			% (len(codon_lst1), len(codon_lst2))
		)
	else:
		length = len(codon_lst1)
	purine = ("A", "G")
	pyrimidine = ("T", "C")
	base_tuple = ("A", "T", "C", "G")
	codon_dict = codon_table.forward_table
	stop = codon_table.stop_codons
	codon_npath = {}
	for i, j in zip(codon_lst1, codon_lst2):
		if i != "---" and j != "---":
			codon_npath.setdefault((i, j), 0)
			codon_npath[(i, j)] += 1
	S_sites = N_sites = 0
	freqSN = [
		{"A": 0, "T": 0, "C": 0, "G": 0},  # synonymous
		{"A": 0, "T": 0, "C": 0, "G": 0},
	]  # nonsynonymous
	for codon_pair, npath in codon_npath.items():
		codon = codon_pair[0]
		if codon not in ['TGA', 'TAA', 'TAG'] and "N" not in codon:
			S = N = 0
			for pos in range(3):
				for base in base_tuple:
					if codon[pos] == base:
						continue
					neighbor_codon = codon[:pos] + base + codon[pos + 1 :]
					if neighbor_codon in stop:
						continue
					weight = pi[neighbor_codon]
					if codon[pos] in pyrimidine and base in pyrimidine:
						weight *= k
					elif codon[pos] in purine and base in purine:
						weight *= k
					if codon_dict[codon] == codon_dict[neighbor_codon]:
						S += weight
						freqSN[0][base] += weight * npath
					else:
						N += weight
						freqSN[1][base] += weight * npath
			S_sites += S * npath
			N_sites += N * npath
	norm_const = 3 * length / (S_sites + N_sites)
	S_sites *= norm_const
	N_sites *= norm_const
	for i in freqSN:
		norm_const = sum(i.values())
		for b in i:
			i[b] /= norm_const
	return S_sites, N_sites, freqSN 

def _get_TV(codon_lst1, codon_lst2, codon_table=default_codon_table): 
	"""Get TV (PRIVATE). 
   
	Arguments: 
	- T - proportions of transitional differences 
	- V - proportions of transversional differences 

	""" 
	purine = ("A", "G") 
	pyrimidine = ("C", "T") 
	TV = [0, 0] 
	sites = 0 
	for codon1, codon2 in zip(codon_lst1, codon_lst2): 
		if "---" not in (codon1, codon2): 
			for i, j in zip(codon1, codon2): 
				if i == j:
					pass 
				elif i in purine and j in purine: 
					TV[0] += 1 
				elif i in pyrimidine and j in pyrimidine: 
					TV[0] += 1 
				else: 
					TV[1] += 1 
				sites += 1 

	return (float(TV[0]) / sites, float(TV[1]) / sites) 

def _get_kappa_t(pi, TV, t=False): 
	"""Calculate kappa (PRIVATE). 
	
	The following formula and variable names are according to PMID: 10666704 
	""" 
	pi["Y"] = pi["T"] + pi["C"] 
	pi["R"] = pi["A"] + pi["G"] 
	A = (
		2 * (pi["T"] * pi["C"] + pi["A"] * pi["G"])
		+ 2
		* (
			pi["T"] * pi["C"] * pi["R"] / pi["Y"]
			+ pi["A"] * pi["G"] * pi["Y"] / pi["R"]
		)
		* (1 - TV[1] / (2 * pi["Y"] * pi["R"]))
		- TV[0]
	) / (2 * (pi["T"] * pi["C"] / pi["Y"] + pi["A"] * pi["G"] / pi["R"]))
	B = 1 - TV[1] / (2 * pi["Y"] * pi["R"]) 
	a = -0.5 * log(A)  # this seems to be an error in YANG's original paper 
	b = -0.5 * log(B) 
	kappaF84 = (a/b) - 1 
	if t is False: 
		kappaHKY85 = 1 + (pi["T"] * pi["C"] / pi["Y"] + pi["A"] * pi["G"] / pi["R"]) * kappaF84 / (pi["T"] * pi["C"] + pi["A"] * pi["G"]) 
		return kappaHKY85 
	else: 
		t = (4 * pi["T"] * pi["C"] * (1 + kappaF84 / pi["Y"]) +  4 * pi["A"] * pi["G"] * (1 + kappaF84 / pi["R"]) + 4 * pi["Y"] * pi["R"]) * b 
		return t
 

def _get_codon_fold(codon_table):
	"""Classify different position in a codon into different folds (PRIVATE)."""

	def find_fold_class(codon, forward_table):
		base = {"A", "T", "C", "G"}
		fold = ""
		codon_base_lst = list(codon)
		for i, b in enumerate(codon_base_lst):
			other_base = base - set(b)
			aa = []
			for j in other_base:
				codon_base_lst[i] = j
				try:
					aa.append(forward_table["".join(codon_base_lst)])
				except KeyError:
					aa.append("stop")
			if aa.count(forward_table[codon]) == 0:
				fold += "0"
			elif aa.count(forward_table[codon]) in (1, 2):
				fold += "2"
			elif aa.count(forward_table[codon]) == 3:
				fold += "4"
			else:
				raise RuntimeError("Unknown Error, cannot assign the position to a fold")
			codon_base_lst[i] = b
		return fold
	fold_table = {}
	for codon in codon_table.forward_table:
		if "U" not in codon:
			fold_table[codon] = find_fold_class(codon, codon_table.forward_table)
	fold_table["---"] = "---"
	return fold_table

def _yn00(seq1, seq2, codon_table=default_codon_table, estimate_kappa=True, median_kappa=5.20):
	"""YN00 method main function (PRIVATE).
	Nomenclature is according to Yang and Nielsen (2000), PMID 10666704.
	"""
	from collections import defaultdict
	from scipy.linalg import expm

	fcodon = [
		{"A": 0, "G": 0, "C": 0, "T": 0},
		{"A": 0, "G": 0, "C": 0, "T": 0},
		{"A": 0, "G": 0, "C": 0, "T": 0},
	]
	codon_fold_dict = _get_codon_fold(codon_table)
	fold0_cnt = defaultdict(int)
	fold4_cnt = defaultdict(int)
	fold_all_cnt = defaultdict(int)
	for codon in seq1 + seq2:
		# count sites at different codon position
		if codon not in ['TGA', 'TAA', 'TAG'] and "N" not in codon:
			if (codon != "---") & (codon != 'TGA'):
				fcodon[0][str(codon[0])] += 1
				fcodon[1][str(codon[1])] += 1
				fcodon[2][str(codon[2])] += 1
			# count sites in different degenerate fold class

			fold_num = codon_fold_dict[codon]

			for i, f in enumerate(fold_num):
				if f == "0":
					fold0_cnt[codon[i]] += 1
				if f == "4":
					fold4_cnt[codon[i]] += 1

	f0_total = sum(fold0_cnt.values())
	f4_total = sum(fold4_cnt.values())

	for i, j in zip(fold0_cnt, fold4_cnt):
		fold0_cnt[i] = fold0_cnt[i] / float(f0_total)
		fold4_cnt[i] = fold4_cnt[i] / float(f4_total)

	if estimate_kappa == True:

		# TODO:
		# the initial kappa is different from what yn00 gives,
		# try to find the problem.
		TV = _get_TV(seq1, seq2, codon_table=codon_table)
		if 0 in TV:
			pass
		else:
			k04 = (_get_kappa_t(fold0_cnt, TV), _get_kappa_t(fold4_cnt, TV))

			kappa = (f0_total * k04[0] + f4_total * k04[1]) / (f0_total + f4_total)
			
			# count synonymous sites and non-synonymous sites
			for i in range(3):
				tot = sum(fcodon[i].values())
				fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
			pi = defaultdict(int)
			for i in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
				if "U" not in i:
					pi[i] = 0
			for i in seq1 + seq2:
				pi[i] += 1
			S_sites1, N_sites1, bfreqSN1 = _count_site_YN00(seq1, seq2, pi, k=kappa, codon_table=codon_table)
			S_sites2, N_sites2, bfreqSN2 = _count_site_YN00(seq2, seq1, pi, k=kappa, codon_table=codon_table)
			N_sites = (N_sites1 + N_sites2) / 2
			S_sites = (S_sites1 + S_sites2) / 2

			return([kappa, S_sites, N_sites])
	if estimate_kappa == False:
		kappa = median_kappa
		# count synonymous sites and non-synonymous sites
		for i in range(3):
			tot = sum(fcodon[i].values())
			fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
		pi = defaultdict(int)
		for i in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
			if "U" not in i:
				pi[i] = 0
		for i in seq1 + seq2:
			pi[i] += 1
		S_sites1, N_sites1, bfreqSN1 = _count_site_YN00(seq1, seq2, pi, k=kappa, codon_table=codon_table)
		S_sites2, N_sites2, bfreqSN2 = _count_site_YN00(seq2, seq1, pi, k=kappa, codon_table=codon_table)
		N_sites = (N_sites1 + N_sites2) / 2
		S_sites = (S_sites1 + S_sites2) / 2		

		return([kappa, S_sites, N_sites])

def pairwise_ts_tv(fasta_dict):

		species = fasta_dict.keys()
		samples = random.sample(species, k=50)

		pairs_strains = combinations(samples,2)

		ratios = list()
		synonymous_sites = list()
		non_synonymous_sites = list()

		for pair in pairs_strains:

			#print(pair)
			seq1 = fasta_dict[pair[0]]
			seq2 = fasta_dict[pair[1]]

			codon1 = [seq1[i:i+3] for i in range(0, len(seq1), 3)]
			codon2 = [seq2[i:i+3] for i in range(0, len(seq2), 3)]

			# Delete stop codons:
			del codon1[-1]
			del codon2[-1]

			kappa = _yn00(codon1, codon2)
			if kappa > 0:
				ratios.append(kappa[0])
				synonymous_sites.append(kappa[1])
				non_synonymous_sites.append(kappa[2])


			# print 'Ts/Tv: %0.3f %d/%d Transitions: A-G:%d C-T:%d Transversions: A-C:%d A-T:%d C-G:%d G-T:%d\n' %(float(Ts)/Tv,Ts,Tv,AG,CT,AC,AT,CG,GT);
		mean_kappa = np.mean(filter(None, ratios))
		mean_synonymous = np.mean(filter(None, synonymous_sites))
		mean_non_synonymous = np.mean(filter(None, non_synonymous_sites))
		print (" This is the mean kappa:")
		print mean_kappa
		print(" this is mean synonymous and non and syn percentage")
		print(mean_synonymous, mean_non_synonymous, mean_synonymous/(mean_synonymous+mean_non_synonymous))
		return([mean_kappa, mean_synonymous, mean_non_synonymous])

def write_table_hist(gene_names, vector, name, results_dir=results_dir):
	d = {"gene_names":gene_names, name: vector}
	print(d)
	df = pd.DataFrame(d)
	df.columns = [name, 'Gene']
	df = df.dropna(thresh=2)
	print(df)

	df.to_csv(results_dir+"{name}_counts.csv".format(name=name))

	print df[name]
	# Histogram of kappa distribution
	f = plt.hist(df[name], bins='auto')  # arguments are passed to np.histogram
	plt.title("Histogram {name}".format(name=name))
	plt.show()

def calling_kappa(filtered_genes_file, group_alignment_dir):
	file_name = filtered_genes_file
	gene_stats = pd.read_table(file_name, delimiter=",")
	candidates = gene_stats["Gene.group"].tolist()

	gene_alignments = [(f) for f in glob.glob(group_alignment_dir+"*.fna")]

	gene_names = list()
	mean_kappa = list()
	mean_synonymous = list()
	mean_nonsynonymous = list()
	counts = 0
	for i in gene_alignments[0:100]:
		base = os.path.basename(i)[:-4]
		if base in candidates: #### I wanna run kappa for every gene 
			counts +=1
			dict_fasta =  parse_fasta(i, basename=base)
			results = pairwise_ts_tv(dict_fasta)
			mean_kappa.append(results[0])
			mean_synonymous.append(results[1])
			mean_nonsynonymous.append(results[2])
			gene_names.append(base)
			print(counts)

	write_table_hist(gene_names, mean_kappa, "Kappa")
	write_table_hist(gene_names, mean_synonymous, "Synonymous_sites")
	write_table_hist(gene_names, mean_nonsynonymous, "Non_Synonymous_sites")

calling_kappa(filtered_genes_file=filtered_genes_directory , group_alignment_dir=group_alignment_dir)


def taking_mean_synonymous_non_synonymous_counts(directory_results=results_dir):

	results = [(f) for f in glob.glob(directory_results+"*.sites")]

	synonynous_counts = 0
	non_synonymous_counts = 0

	#sitesS	sitesN
	for i in results:
		df = pd.read_csv(i, delimiter=",")
		synonynous_counts += df["sitesS"].mean()
		non_synonymous_counts += df["sitesN"].mean()

	print(synonynous_counts, non_synonymous_counts)

#taking_mean_synonymous_non_synonymous_counts()