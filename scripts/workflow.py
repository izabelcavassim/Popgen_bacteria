# Copyright 2020 by Maria Izabel Cavassim Alves (izabelcavassim@gmail.com).
# All rights reserved.
# This code is part of the workflow for the manuscript "Recombination facilitates adaptive evolution in rhizobial soil bacteria"
# It uses the gwf workflow software see: https://gwf.app/
# It is written in python 3. 
# Scripts used in this pipeline our found in: https://github.com/izabelcavassim/Popgen_bacteria/tree/master/scripts
# Directories are merely representative and should be changed if scripts are used
# Feel free to contact me with questions!

from gwf import Workflow
import glob
import os
import itertools
gwf = Workflow()

# running clonal-frame

def spliting_fasta_files():
	outputs=['/home/mica16/NChain/faststorage/rhizobium/SFS_data/gsA/group1000.fna', '/home/mica16/NChain/faststorage/rhizobium/SFS_data/gsB/group1000.fna', '/home/mica16/NChain/faststorage/rhizobium/SFS_data/gsC/group1000.fna', '/home/mica16/NChain/faststorage/rhizobium/SFS_data/gsD/group1000.fna', '/home/mica16/NChain/faststorage/rhizobium/SFS_data/gsE/group1000.fna']
	options = {
	'memory':'8g',
	'cores':'8',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = '''
	python clonal_frame_analysis.py
	'''
	return inputs, outputs, options, spec

def Clonal_frame(gene_name, genospecies):
	inputs = []
	outputs = [f'{gene_name}_{geno}_cf_final.tab']
	options = {
	'memory':'1gb',
	'cores':'1',
	'walltime':'20:00',
	'account': 'NChain'
	}
	spec = f'''
	ClonalFrameML tree_{genospecies}.nw /home/mica16/NChain/faststorage/rhizobium/SFS_data/{genospecies}/{gene_name} {gene_name}_{genospecies}_clonal_frame -kappa 5.20 -emsim 1000  > {gene_name}_clonal_frame_logs.txt
	grep -E "(^nu)|(^1/delta)|(^R/theta)" {gene_name}_{genospecies}_clonal_frame.em.txt | awk '{{ print "{gene_name}\t{genospecies}\t" $0}}' > {gene_name}_{geno}_cf_final.tab
		
	mv {gene_name}_{genospecies}_clonal_frame.emsim.txt output_emsim/
	mv {gene_name}_{genospecies}_clonal_frame.importation_status.txt output_importation/ 
	mv {gene_name}_{genospecies}_clonal_frame.ML_sequence.fasta output_ML_fasta/ 
	mv {gene_name}_{genospecies}_clonal_frame.labelled_tree.newick output_ML_fasta/ 

	cat {gene_name}_{geno}_cf_final.tab >> final_clframe.tab

	'''
	return inputs, outputs, options, spec

def transition_transversion_analysis(script, specific_genes = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", 
											 gene_alignments_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/group_alns/",
											 results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/"):
	inputs = [f'{script}', "/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv"]
	outputs = [f'{results_dir}'+"Non_Synonymous_sites_counts.csv", f'{results_dir}'+"Synonymous_sites_counts.csv"]
	options = {
	'memory':'1gb',
	'cores':'1',
	'walltime':'12:00:00',
	'account': 'NChain'
	}
	spec = f'''python ts_tv_calc.py {specific_genes} {gene_alignments_dir} {results_dir}'''
	print(spec)
	return inputs, outputs, options, spec

def SFS_per_geno(script, filtered_name_syn ="/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", 
											 metadata = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Rhizobium_soiltypes_new.txt",
											 results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/",
											 snp_file = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/newsnps_100.hdf5"):
	inputs = ["/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv"]
	outputs = [f'{results_dir}'+"SFS_gsE_counts.csv"]
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = f'''
		source activate py36
		python SFS_per_geno.py {filtered_name_syn} {metadata} {results_dir} {snp_file}'''
	print(spec)
	return inputs, outputs, options, spec


def LD_per_geno_per_gene(script, filtered_name_syn = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", 
											 metadata = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Rhizobium_soiltypes_new.txt",
											 results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/",
											 snp_file = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/newsnps_100.hdf5"):
	inputs = ["/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv"]
	outputs = [f'{results_dir}LD_decay_gsE_average_r2_10_snps_min_maf_0.01.csv']
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = f'''
		source activate py36
		python LD_per_gene.py {filtered_name_syn} {metadata} {results_dir} {snp_file}'''
	print(spec)
	return inputs, outputs, options, spec


def LD_per_geno(script, metadata = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Rhizobium_soiltypes_new.txt",
				results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/",
				snp_file = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/newsnps_100.hdf5",
				filtered_genes = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", 
				maf = 0.01, bin_size = 10):
	inputs = [ "/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv"]
	outputs = [f'{results_dir}plotting_intergenic_LD_0_{bin_size}_gsE.csv']
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = f'''
		source activate py36
		python Intragenic_LD.py {filtered_genes} {metadata} {results_dir} {snp_file} {maf} {bin_size}'''
	print(spec)
	return inputs, outputs, options, spec

def Grapes_preparation(script, filtered_genes="/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", filtered_name_syn = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Synonymous_sites_counts.csv",
											 filtered_name_non_syn = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Non_Synonymous_sites_counts.csv",
											 metadata = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Rhizobium_soiltypes_new.txt",
											 results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/",
											 snp_file = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/newsnps_100.hdf5", quantiles=False,input_recombination_files="/home/mica16/NChain/faststorage/rhizobium/SFS_data/grapes_files/"):
	inputs = ["/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Non_Synonymous_sites_counts.csv", "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Synonymous_sites_counts.csv"]
	genos = ['gsA', 'gsB','gsC', 'gsD', 'gsE']
	comb = list(itertools.permutations(genos, 2))
	out = list()
	for c in comb:
		out.append(f'{results_dir}'+f'{c[0]}_{c[1]}_grapes.txt')
	outputs = out
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = f'''
		source activate py36
		python Grapes_files_pairwise_quantile.py {filtered_genes} {filtered_name_syn} {filtered_name_non_syn} {metadata} {results_dir} {snp_file} {quantiles} {input_recombination_files}'''
	print(spec)
	return inputs, outputs, options, spec


def Grapes_preparation_quantiles(script, filtered_genes="/home/mica16/NChain/faststorage/rhizobium/SFS_data/filtered_gene_statistics.csv", filtered_name_syn = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Synonymous_sites_counts.csv",
											 filtered_name_non_syn = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Non_Synonymous_sites_counts.csv",
											 metadata = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Rhizobium_soiltypes_new.txt",
											 results_dir =  "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/",
											 results_dir_2 = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/",
											 snp_file = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/newsnps_100.hdf5", quantiles=False, input_recombination_files="/home/mica16/NChain/faststorage/rhizobium/SFS_data/grapes_files/"):
	inputs = ["/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Non_Synonymous_sites_counts.csv", "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Synonymous_sites_counts.csv"]
	genos = ['gsA', 'gsB','gsC', 'gsD', 'gsE']
	comb = list(itertools.permutations(genos, 2))
	out = list()
	qs = ["first", "second", "third"]
	for c in comb:
		for q in qs:
			out.append(f'{results_dir_2}'+f'{c[0]}_{c[1]}__{q}_quantile_grapes.txt')
	outputs = out
	#print(outputs)
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'00:20:00',
	'account': 'NChain'
	}
	spec = f'''
		source activate py36
		python Grapes_files_pairwise_quantile.py {filtered_genes} {filtered_name_syn} {filtered_name_non_syn} {metadata} {results_dir} {snp_file} {quantiles} {input_recombination_files}'''
	print(spec)
	return inputs, outputs, options, spec

def Grapes_run_full(results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/", geno1=None, geno2=None):
	inputs = [f'{results_dir}'+f'{geno1}_{geno2}_grapes.txt']
	outputs = [f'{results_dir}Grapes_results/{geno1}_{geno2}_grapes_output.txt']
	spec = f'''
	source activate py36
	/home/mica16/grapes-izabel-work/bin/grapes -in {results_dir}{geno1}_{geno2}_grapes.txt -out {results_dir}Grapes_results/{geno1}_{geno2}_grapes_output.txt -nb_rand_start 20 -model all
	'''
	options = {
	'memory':'8gb',
	'cores':'1',
	'walltime':'24:00:00',
	'account': 'NChain'
	}

	return inputs, outputs, options, spec

def Grapes_run_quartiles(quantile, results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Grapes_results/Quantiles/", inputs_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/", geno1=None, geno2=None):
	#gsB_gsE_first_quantile_grapes_output_quantiles.txt
	inputs = [f'{inputs_dir}{geno1}_{geno2}__{quantile}_quantile_grapes.txt']
	outputs = [f'{results_dir}{geno1}_{geno2}_{quantile}_quantile_grapes_output_quantiles.txt']
	spec = f'''
	source activate /home/mica16/anaconda2/envs/py36
	/home/mica16/grapes-izabel-work/bin/grapes -in {inputs_dir}{geno1}_{geno2}__{quantile}_quantile_grapes.txt -out {results_dir}{geno1}_{geno2}_{quantile}_quantile_grapes_output_quantiles.txt -nb_rand_start 20 -model all
	'''
	options = {
	'memory':'20gb',
	'cores':'8',
	'walltime':'20:00:00',
	'account': 'NChain'
	}

	return inputs, outputs, options, spec


def Grapes_run_quartiles_simulations(quantile_file, name, results_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Grapes_results/Quantiles_simulations/"):
	#gsB_gsE_first_quantile_grapes_output_quantiles.txt
	inputs = [quantile_file]
	outputs = [f'{results_dir}out_{name}']
	spec = f'''
	source activate /home/mica16/anaconda2/envs/py36
	/home/mica16/grapes-izabel-work/bin/grapes -in {quantile_file} -out {results_dir}out_{name} -model GammaZero
	'''
	options = {
	'memory':'20gb',
	'cores':'8',
	'walltime':'20:00:00',
	'account': 'NChain'
	}

	return inputs, outputs, options, spec

def Run_enfcprime(gene_name, gene_directory, output_dir, n, geno):
	inputs = [f'{gene_directory}']
	outputs = [f'{output_dir}{geno}_{gene_name}']
	spec = f'''
	/home/mica16/NChain/faststorage/rhizobium/SFS_data/ENCprime/bin/SeqCount -c {gene_directory} {n}
	/home/mica16/NChain/faststorage/rhizobium/SFS_data/ENCprime/bin/ENCprime {gene_directory}.codcnt /home/mica16/NChain/faststorage/rhizobium/SFS_data/{geno}/concatenated_{geno}.fasta.acgtfreq 11 {geno}_{gene_name} 1 -q
	mv {geno}_{gene_name} {output_dir}
	'''
	options = {
	'memory':'1gb',
	'cores':'8',
	'walltime':'20:00:00',
	'account': 'NChain'
	}

	return inputs, outputs, options, spec

############################################################################
# Running analyses by calling the functions above						   #
############################################################################

genospecies = ["gsA", "gsB", "gsC", "gsD", "gsE"]

# Running clonal frame per gene
for geno in genospecies:
 	genes = [(f) for f in glob.glob("/home/mica16/NChain/faststorage/rhizobium/SFS_data/{geno}/*.fna".format(geno=geno))]
 	for gene in genes:
 		gene = os.path.basename(gene)
 		gwf.target_from_template("Clonal_frame_{geno}_{gene}".format(geno=geno, gene=gene), Clonal_frame(gene_name=gene, genospecies=geno))

gwf.target_from_template("Transition_transversion_all", transition_transversion_analysis(script="ts_tv_calc.py"))
gwf.target_from_template("SFS_each_geno", SFS_per_geno(script="SFS_per_geno.py"))
gwf.target_from_template("LD_each_gene", LD_per_geno(script="Intragenic_LD.py"))
gwf.target_from_template("LD_each_geno", LD_per_geno_per_gene(script="Intragenic_LD.py"))
gwf.target_from_template("Grapes_input_files_pairwise", Grapes_preparation(script="Grapes_files_pairwise_quantile.py", input_recombination_files="/home/mica16/NChain/faststorage/rhizobium/SFS_data/grapes_files/"))
gwf.target_from_template("Grapes_input_files_quantiles", Grapes_preparation_quantiles(script="Grapes_files_pairwise_quantile.py", quantiles=True,  results_dir="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/",  results_dir_2="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/", input_recombination_files="/home/mica16/NChain/faststorage/rhizobium/SFS_data/grapes_files/"))
gwf.target_from_template("Grapes_input_files_quantiles_ClonalFrame", Grapes_preparation_quantiles(script="Grapes_files_pairwise_quantile.py", quantiles=True, results_dir="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination_ClonalFrame/",  results_dir_2="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination_ClonalFrame/", input_recombination_files="/home/mica16/NChain/faststorage/rhizobium/SFS_data/grapes_files_ClonalFrame/"))

# Calling grapes pairwise estimates of alpha
genos = ['gsA', 'gsB','gsC', 'gsD', 'gsE']
comb = list(itertools.permutations(genos, 2))
for c in comb:
	gwf.target_from_template(f'Grapes_run_{c[0]}_{c[1]}', Grapes_run_full(geno1=c[0], geno2=c[1]))

# Calling grapes for different recombination quantiles
genos = ['gsA','gsB','gsC','gsD','gsE']
quantiles = ["first", "second", "third"] #, "fourth"]
for c in comb:
	for q in quantiles:
		gwf.target_from_template(f'Grapes_run_{c[0]}_{c[1]}_quantile_{q}', Grapes_run_quartiles(geno1=c[0], geno2=c[1], quantile=q, inputs_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/", results_dir="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Grapes_results/Quantiles/"))
		gwf.target_from_template(f'Grapes_run_{c[0]}_{c[1]}_quantile_{q}_ClonalFrame', Grapes_run_quartiles(geno1=c[0], geno2=c[1], quantile=q, inputs_dir = "/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination_ClonalFrame/", results_dir="/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Grapes_results/Quantiles_ClonalFrame/"))


# Running simulations by shuffling recombination classes 
grapes_quantiles_files = [(f) for f in glob.glob("/home/mica16/NChain/faststorage/rhizobium/SFS_data/Results/Quartiles_files_recombination/*.txt")]
for file in grapes_quantiles_files:
	name = os.path.basename(file)
	gwf.target_from_template(f'{name}', Grapes_run_quartiles_simulations(quantile_file=file, name=name))