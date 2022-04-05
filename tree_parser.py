"""

author: Michael Dubovsky, Technion Israel Institute of Technology @

moodule to parse IgPhyML phylogenetic tree output and store the information output from Tree Analyzer Block.

"""
################################################################################
#Imports Section
import os
import tree_analyzer as tree_functions
import time 
import pandas as pd
import matplotlib.pyplot as plt
import argparse
################################################################################
#Declare Global Variables
start_time = time.time()
df_results = pd.DataFrame()
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
#Spike Indicator Decision for Clone Lineage
#Function to decide the spike clone indicator by crossing with metadata from AIRR DB file from ImmuneDB.
#"+" = Spike+; "-" = Spike-; "+/-" = Mixed Spike+ & Spike-
def clone_spike_indicator(sequence_id,df_immune_db):
	spike_plus_indicator = 0
	spike_minus_indicator = 0
	for id in sequence_id:
		#####################################################################################
		#Handle the issue that IgPhyML flow changed the sequence ID character ":" to "-" (issue with IgPhyML flow!)
		id_original = id.replace("-",":")
		id_original = id_original.replace("M03592:154:000000000:","M03592:154:000000000-")
		#####################################################################################
		airr_row_seq_id = df_immune_db.loc[df_immune_db['sequence_id'] == id_original]
		if "Spike+" in airr_row_seq_id.at[int(airr_row_seq_id.index.values),'METADATA_cell_subset_description']:
			spike_plus_indicator = 1
		else:
			spike_minus_indicator = 1
	if (spike_minus_indicator) and (spike_plus_indicator):
		return "+/-"
	elif (spike_minus_indicator):
		return "-"
	elif (spike_plus_indicator):
		return "+"
	else :
		return ""
################################################################################
#Parse Command Line Argument Section
#Create object for parsing command-line options
parser = argparse.ArgumentParser(description="Command line for Tree Parser,Reportoire file,Fasta Files Directory")
parser.add_argument("-r", "--repertoire", type=str, help="Path to IgPhyML repertoire file")
parser.add_argument("-f", "--fasta", type=str, help="Path to IgPhyML FASTA files DIR")
parser.add_argument("-o", "--output", type=str, help="Output CSV name")
parser.add_argument("-a", "--airr_db_file", type=str, help="Path to AIRR DB file for Meta Data extraction")
#Parse the command line arguments to an object
args = parser.parse_args()
#Safety if no parameter have been given
if not args.repertoire:
	print("No Repertoire File has been given. Please insert reperotire file in command line --repertoire")
	print("For help type --help")
	exit()
if not args.fasta:
	print("No FASTA Files DIR has been given. Please insert FASTA File DIR in command line --fasta")
	print("For help type --help")
	exit()
if not args.output:
	print("No Output File Name has been given. Please insert Output File Name in command line --output")
	print("For help type --help")
	exit()
if not args.airr_db_file:
	print("No AIRR DB File has been given. Please insert AIRR DB File in command line --airr_db_file")
	print("For help type --help")
	exit()
################################################################################
#Inject original Immune DB table to Panda
df_immune_db = pd.read_csv(args.airr_db_file, sep='\t')
################################################################################
#Read Repertoire Trees File from IgphyML
repertoire_file = open(args.repertoire, "r")
################################################################################
#Parse Repertoire Trees File from IgphyML
start_parsing_flag = 0
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
	if (start_parsing_flag > 2):
                split_line = line.split()
                clone_number = split_line[0]
                tree_topology_str = split_line[14]
                fasta_file_str = str(clone_number) + ".fasta"
                with open((args.fasta + "/" + fasta_file_str) , "r") as fasta_file:
                        fasta_str = fasta_file.read()
                        fasta_file.close()
                clone_results_container = tree_functions.tree_analyzer_flow(fasta_str, tree_topology_str, df_immune_db)
                clone_results_container['clone_id'] = clone_number
                clone_results_container['spike']= clone_spike_indicator(clone_results_container['sequence_id_list'], df_immune_db)
		#Move results to Panda Table !!!
                df_results = df_results.append(clone_results_container, ignore_index=True)
repertoire_file.close()
##Dump the data to final csv file
df_results.to_csv(args.output, index=False)
################################################################################
print("--- %s seconds ---" % (time.time() - start_time))
################################################################################
