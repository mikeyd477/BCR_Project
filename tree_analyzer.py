# -*- coding: utf-8 -*-
"""

"""

################################################################################
import os
#os.system("pip install ete3 pyqt5")
#os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree
from ete3 import PhyloTree
#from ete3 import TreeStyle
from collections import Counter
################################################################################
#Definitions
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
#Distance in Mutation calculation function.
#Function to calculate distance in mutation between two given neclutide sequences.
def mutation_dist_calc(germ_seq,seq):
  distance = 0
  # Iterate over the sequences
  if (len(germ_seq) != len(seq)):
    print("Error the alignment sequence and the germline seq are with different sizes!!!")
  for index in range(0, len(germ_seq)):
    if (germ_seq[index] == "N" ):
      continue
    if (germ_seq[index] == "-"):
      continue
    if (seq[index] == "-"):
      continue
    if (seq[index] == "N" ):
      continue
    if (germ_seq[index] != seq[index]):
      distance = distance + 1
  return distance
################################################################################
#CDR3 Consensus Handling
#Function to insert the CDR3 consensus by checking the majority between the given sequences.
def cdr3_consensus_calc(germ_seq,seqs):
  germline_sequence = ""
  for index in range(0, len(germ_seq)):
    if (germ_seq[index] != "-"):
      germline_sequence = germline_sequence + germ_seq[index]
      continue
    else:
      neuc_string = ""
      for neuc in seqs:
        neuc_string = neuc_string + neuc[index]
      consensus = Counter(neuc_string)
      consensus = max(consensus, key = consensus.get)
      germline_sequence = germline_sequence + str(consensus)
  return germline_sequence
################################################################################
#Synonymous non Synonymous mutations counter
#Function to count synonymous non synonymous between germline sequence and SUT, sequence under test.
#Please note that in this function we calculate synonymous and not synonymous mutations between
#Amino Acids : V,A,T,G,P (the AA that have 4 codons coding to them).
def synonymous_mutation_calc(germ_seq,seq):
  count_synon = 0
  count_not_synon = 0
  start_calc_flag = 0
  max_synonymous_mutations = 0
  for sequence in seq:
    count_index = 1
    count_synon_mutation_in_sequence = 0
    codon_str_seq = ""
    codon_str_germ = ""
    aa_str_seq = ""
    aa_str_germ = ""
    for index in range(0, len(sequence)):
      if (sequence[index] != "N"):
        start_calc_flag = 1
      if (start_calc_flag == 0):
        continue
      codon_str_seq = codon_str_seq + sequence[index]
      codon_str_germ = codon_str_germ + germ_seq[index]
      if ((count_index % 3) == 0):
        #one of the neclutides is N therefore not possible to deduce
        if (("N" in codon_str_seq) or ("N" in codon_str_germ)): 
          codon_str_seq = ""
          codon_str_germ = ""
          count_index = 1
          continue
        if (("-" in codon_str_seq) or ("-" in codon_str_germ)): 
          codon_str_seq = ""
          codon_str_germ = ""
          count_index = 1
          continue
        #No mutations in codon
        if (codon_str_seq == codon_str_germ):
          codon_str_seq = ""
          codon_str_germ = ""
          count_index = 1
          continue
        #Translate codons to AA    
        aa_str_seq = aa_dict[codon_str_seq]
        aa_str_germ = aa_dict[codon_str_germ]
        if ((aa_str_germ == "V") or (aa_str_germ == "A") or (aa_str_germ == "T") or (aa_str_germ == "G") or (aa_str_germ == "P")):
          if  (aa_str_germ == aa_str_seq):
            count_synon = count_synon + 1
            count_synon_mutation_in_sequence = count_synon_mutation_in_sequence +1
          else:
            count_not_synon = count_not_synon + 1
        codon_str_seq = ""
        codon_str_germ = ""
        count_index = 1
      else:
        count_index = count_index + 1
    if (count_synon_mutation_in_sequence > max_synonymous_mutations):
        max_synonymous_mutations = count_synon_mutation_in_sequence
  return count_synon,count_not_synon,max_synonymous_mutations
################################################################################
#Synonymous non Synonymous mutations counter
#Function to count synonymous non synonymous mutation ratio between germline sequence and farthest sequence and closest sequence to germline.
#Please note that in this function we calculate synonymous and not synonymous mutations for all amino acids.
def synonymous_mutation_ratio_per_seq_calc(germ_seq,sequence):
  count_synon = 0
  count_not_synon = 0
  start_calc_flag = 0
  count_index = 1
  codon_str_seq = ""
  codon_str_germ = ""
  aa_str_seq = ""
  aa_str_germ = ""
  for index in range(0, len(sequence)):
    if (sequence[index] != "N"):
        start_calc_flag = 1
    if (start_calc_flag == 0):
        continue
    codon_str_seq = codon_str_seq + sequence[index]
    codon_str_germ = codon_str_germ + germ_seq[index]
    if ((count_index % 3) == 0):
    #one of the neclutides is N therefore not possible to deduce
        if (("N" in codon_str_seq) or ("N" in codon_str_germ)): 
            codon_str_seq = ""
            codon_str_germ = ""
            count_index = 1
            continue
        if (("-" in codon_str_seq) or ("-" in codon_str_germ)): 
            codon_str_seq = ""
            codon_str_germ = ""
            count_index = 1
            continue
        #No mutations in codon
        if (codon_str_seq == codon_str_germ):
            codon_str_seq = ""
            codon_str_germ = ""
            count_index = 1
            continue
        #Translate codons to AA    
        aa_str_seq = aa_dict[codon_str_seq]
        aa_str_germ = aa_dict[codon_str_germ]
        if  (aa_str_germ == aa_str_seq):
            count_synon = count_synon + 1
        else:
            count_not_synon = count_not_synon + 1
        codon_str_seq = ""
        codon_str_germ = ""
        count_index = 1
    else:
        count_index = count_index + 1
  if (count_synon == 0):
    sequence_omega = 0
  else:
    sequence_omega = (count_not_synon / count_synon)
  return sequence_omega
################################################################################
#Synonymous non Synonymous mutations counter
#Function to count synonymous non synonymous between germline sequence and SUT, sequence under test.
def tree_analyzer_flow(fasta,tree_topology_str):
    # Load a tree and link it to an alignment.
    tree_topolo_str = "((((((M03592-154-000000000-JCY9V-1-2114-26424-16841:0.0558811845,(M03592-154-000000000-JCY9V-1-2108-9799-13897:0.0104575587,M03592-154-000000000-JCY9V-1-2117-9143-9891:0.0000000100):0.0209073238):0.0000008340,M03592-154-000000000-JCY9V-1-2105-2262-11199:0.0000000100):0.0313860135,M03592-154-000000000-JCY9V-1-1115-19386-20227:0.0104092512):0.0527013194,M03592-154-000000000-JCY9V-1-2108-23373-15593:0.0670429058):0.0103785606,M03592-154-000000000-JCY9V-1-2116-13304-13947:0.0000000100):0.1229925141,615860_GERM:0.0000100000);"
    t = PhyloTree(tree_topology_str, alignment=fasta, alg_format="fasta")
	#list to hold all sequence ID's
    sequence_id_list = []
################################################################################
    one_sequence_flag = 0
    #Check if we have 1 actual sequence in the clone to handle the sequence multiplication in igphyml flow.
    for node in t.traverse(strategy="postorder"):
        if node.is_leaf():
            name_len = len(node.name)
            if (node.name[name_len-2:] == "_1"):
                one_sequence_flag = 1
                break
################################################################################
    #Count number of leafs by checking if the branch distance is higher than the threshold.
    #Insert Leafs sequences to List
    #Currently we've set thershold of e-4. that is ~0.03 mutations in our case. 
    count_leafs = 0
    if (one_sequence_flag == 1):
        count_leafs = 1
    else:
        leaf_nodes_list = []
        for node in t.traverse(strategy="postorder"):
            if node.is_leaf():
                if node.dist > 0.0001:
                    count_leafs = count_leafs+1
                    leaf_nodes_list.append(node)
################################################################################
    #Find the germline sequence and insert all other sequences to list:
    germline_sequence = ""
    germline_sequence_with_cdr3 = ""
    sequences_list = []
    for node in t.traverse(strategy="postorder"):
        if node.is_leaf():
            if "GERM" in node.name:
                #print(node)
                #print(node.sequence)
                germline_sequence = node.sequence
            else:
                name_len = len(node.name)
                if (node.name[name_len-2:] != "_1"):
                    sequences_list.append(node.sequence)
                    sequence_id_list.append(node.name)
################################################################################
    #Handle the germline CDR3 consensus:
    germline_sequence_with_cdr3 = cdr3_consensus_calc(germline_sequence,sequences_list)
    #print("Germline Sequence with CDR3 Consensus is")
    #print(germline_sequence_with_cdr3)
    ################################################################################
    #Count the Tree depth by calculating the distance in mutations from the germline!
    #Assumptions:
    #If we have N's(Unkonws) then we assume equal to germline
    #If we have gap(---) then we assume equal to germline
    #Cost function of every neclutide mutation is 1
    max_distance_mutations = 0
    max_distance_seq_name = ""
    max_distance_sequence = ""
    min_distance_mutations = 500
    min_distance_seq_name = ""
    min_distance_sequence = ""
    for node in t.traverse(strategy="postorder"):
        if node.is_leaf():
            if "GERM" in node.name:
                continue
            mutation_distance = mutation_dist_calc(germline_sequence_with_cdr3,node.sequence)
            #print("sequence name ",node.name)
            #print("sequence distance ",mutation_distance)
            if (mutation_distance > max_distance_mutations):
                max_distance_mutations = mutation_distance
                max_distance_seq_name = node.name
                max_distance_sequence = node.sequence
            if (mutation_distance < min_distance_mutations):
                min_distance_mutations = mutation_distance
                min_distance_seq_name = node.name
                min_distance_sequence = node.sequence
################################################################################
    #Count the Max/Min Leaf sequneces distance
    #Assumptions:
    #If we have N's(Unkonws) then we assume equal to germline
    #If we have gap(---) then we assume equal to germline
    #Cost function of every neclutide mutation is 1
    max_distance_mutations_leafs = 0
    min_distance_mutations_leafs = 500
    max_distance_seq_name_leafs_a = ""
    min_distance_seq_name_leafs_a = ""
    max_distance_seq_name_leafs_b = ""
    min_distance_seq_name_leafs_b = ""
    if (one_sequence_flag == 1):
        max_distance_mutations_leafs = 0
        max_distance_seq_name_leafs_a = "NULL"
        max_distance_seq_name_leafs_b = "NULL"
        min_distance_mutations_leafs = 0
        min_distance_seq_name_leafs_a = "NULL"
        min_distance_seq_name_leafs_b = "NULL"
    else:
        for nodeA in leaf_nodes_list:
            for nodeB in leaf_nodes_list:
                if (nodeA.name == nodeB.name):
                    continue
                mutation_distance = mutation_dist_calc(nodeA.sequence,nodeB.sequence)
                if (mutation_distance > max_distance_mutations_leafs):
                    max_distance_mutations_leafs = mutation_distance
                    max_distance_seq_name_leafs_a = nodeA.name
                    max_distance_seq_name_leafs_b = nodeB.name
                if (mutation_distance < min_distance_mutations_leafs):
                    min_distance_mutations_leafs = mutation_distance
                    min_distance_seq_name_leafs_a = nodeA.name
                    min_distance_seq_name_leafs_b = nodeB.name
################################################################################
    #Count the synonymous & non synonymous mutations.
    #Currently counted for AA: V,A,T,G,P
    #Assumptions:
    count_synonymous = 0;
    count_not_synonymous = 0;
    count_synonymous,count_not_synonymous,max_synonymous_mutations_in_sequence = synonymous_mutation_calc(germline_sequence_with_cdr3,sequences_list)
################################################################################
    #Calculate the non synonymous/synonymous mutations ratio for the farthest and closest sequence to germline.
    omega_farthest = synonymous_mutation_ratio_per_seq_calc(germline_sequence_with_cdr3, max_distance_sequence)
    omega_closest = synonymous_mutation_ratio_per_seq_calc(germline_sequence_with_cdr3, min_distance_sequence)
################################################################################
   #Print Results:
    #print("---------Tree Analyzer Results---------")
    #print("Number of Leafs = ",count_leafs)
    #print("Max Mutations Distance Any Sequence ID is ",max_distance_seq_name)
    #print("Max Mutations Distance Any Sequence from Germline is ",max_distance_mutations)
    #print("Min Mutations Distance Any Sequence ID is ",min_distance_seq_name)
    #print("Min Mutations Distance Any Sequence from Germline is ",min_distance_mutations)
    #print("Max Mutations Distance BTW leaves is seqA ID ",max_distance_seq_name_leafs_a, " seqB ID ",max_distance_seq_name_leafs_b)
    #print("Max Mutations Distance BTW leaves is ",max_distance_mutations_leafs)
    #print("Min Mutations Distance BTW leaves is seqA ID ",min_distance_seq_name_leafs_a, " seqB ID ",min_distance_seq_name_leafs_b)
    #print("Min Mutations Distance BTW leaves is ",min_distance_mutations_leafs)
    #print("Number of synonymous mutations is ",count_synonymous)
    #print("Number of not synonymous mutations is ",count_not_synonymous)
    #print("Max Number of synonymous mutations in sequence is ",max_synonymous_mutations_in_sequence)
    #print("Farthest sequence omega ratio is ", omega_farthest)
    #print("Closest sequence omega ratio is ", omega_closest)
################################################################################
    #Insert results to list
    results = {
        "num_leaves": count_leafs,
        "max_mutations_distance_sequence": max_distance_mutations,
        "min_mutations_distance_sequence": min_distance_mutations,
        "max_mutations_between_leaves": max_distance_mutations_leafs,
        "min_mutations_between_leaves": min_distance_mutations_leafs,
        "num_synonymous_mutations_four_fold": count_synonymous,
        "num_not_synonymous_mutations_four_fold": count_not_synonymous,
        "max_synonymous_mutations_in_sequence_four_fold": max_synonymous_mutations_in_sequence,
        "omega_ratio_farthest_sequence": omega_farthest,
        "omega_ratio_closest_sequence": omega_closest,
        "sequence_id_list": sequence_id_list
    }
    return results