################################################################################
import os
import tree_analyzer as tree_functions
import time 
import pandas as pd
################################################################################
#Declare Global Variables
start_time = time.time()
df_results = pd.DataFrame()
################################################################################
#Inject original Immune DB table to Panda
df_immune_db = pd.read_csv("all_samples_together_manipulated_removed_problem_clones.tsv", sep='\t')
print(df_immune_db.head())
################################################################################
#Define Codon 2 Amino Acid Dict
aa_dict = {
  "TTT" : "F",
  "TTC" : "F",
  "TTA" : "L",
  "TTG" : "L",
  "TCT" : "S",
  "TCC" : "S",
  "TCA" : "S",
  "TCG" : "S",
  "TAT"  : "Y",
  "TAC"  : "Y",
  "TAA"  : "Stop",
  "TAG"  : "Stop",
  "TGT"  : "C",
  "TGC"  : "C",
  "TGA"  : "Stop",
  "TGG"  : "W",
  "CTT" : "L",
  "CTC" : "L",
  "CTA" : "L",
  "CTG" : "L",
  "CCT" : "P",
  "CCC" : "P",
  "CCA" : "P",
  "CCG" : "P",
  "CAT" : "H",
  "CAC" : "H",
  "CAA" : "Q",
  "CAG" : "Q",
  "CGT" : "R",
  "CGC" : "R",
  "CGA" : "R",
  "CGG" : "R",
  "ATT" : "I",
  "ATC" : "I",
  "ATA" : "I",
  "ATG" : "M",
  "ACT" : "T",
  "ACC" : "T",
  "ACA" : "T",
  "ACG" : "T",
  "AAT" : "N",
  "AAC" : "N",
  "AAA" : "K",
  "AAG" : "K",
  "AGT" : "S",
  "AGC" : "S",
  "AGA" : "R",
  "AGG" : "R",
  "GTT" : "V",
  "GTC" : "V",
  "GTA" : "V",
  "GTG" : "V",
  "GCT" : "A",
  "GCC" : "A",
  "GCA" : "A",
  "GCG" : "A",
  "GAT" : "D",
  "GAC" : "D",
  "GAA" : "E",
  "GAG" : "E",
  "GGT" : "G",
  "GGC" : "G",
  "GGA" : "G",
  "GGG" : "G"
}
################################################################################
#Read Repertoire Trees File from IgphyML
repertoire_file = open("subject005_sequence_clones_time_point4_igphyml-pass_test.tab", "r")
#with open("615860_all_time_points.fasta", "r") as fasta_file:
################################################################################
#Parse Repertoire Trees File from IgphyML
start_parsing_flag = 0
#print("-------------------------- Test Parsing ----------------------")
while True:
	fasta_str = ""
	fasta_file_str = ""
	split_line = ""
	clone_number = 0
	tree_topology_str = ""
	start_parsing_flag += 1
	line = repertoire_file.readline()
	if not line:
		break
	#print("----------- Line Number ------------- ", start_parsing_flag," -------------------")
	#print(line)
	if (start_parsing_flag > 2):
		split_line = line.split()
		clone_number = split_line[0]
		tree_topology_str = split_line[14]
		#print("---------------------- Test ---------------")
		#print("Clone number is: ", clone_number)
		#print("Tree topology string is: ",tree_topology_str)
		fasta_file_str = str(clone_number) + ".fasta"
		#print("Fasta File string is: ", fasta_file_str)
		with open(("subject005_sequence_clones_time_point4/" + fasta_file_str) , "r") as fasta_file:
			fasta_str = fasta_file.read()
			fasta_file.close()
			#print(fasta_str)
		clone_results_container = tree_functions.tree_analyzer_flow(fasta_str, tree_topology_str)
		clone_results_container['clone_id'] = clone_number
		#print("--------------Test Results Container---------------------")
		#print(clone_results_container)
		################################################################################
		#Move results to Panda Table !!!
		df_results = df_results.append(clone_results_container, ignore_index=True)
		print(df_results.head(6))
		print(df_results[df_results["clone_id"] == "605501"])
repertoire_file.close()
################################################################################
print("--- %s seconds ---" % (time.time() - start_time))
################################################################################
#test_str = file.readline()
#print(test_str)
#test_split = test_str.split()
#print("-------------Test Parsing -------------")
##clone number
#print(test_split[0])
#Tree Topology
#print(test_split[14])
#tree_functions.tree_analyzer_flow(fasta_data, test_split[14])
#file.close()