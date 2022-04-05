"""

author: Michael Dubovsky, Technion Israel Institute of Technology

moodule to analyze features of phylogenetic trees for time dependency study purposes.

The moodule runs with the wrapper "tree_parser"

"""
################################################################################
#Import Section
import os
from ete3 import Tree
from ete3 import PhyloTree
from collections import Counter
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
#Distance in Mutation calculation function.
#Function to calculate distance in mutation between two given neclutide sequences.
def mutation_dist_calc(germ_seq,seq,seq_name):
  distance = 0
  # Iterate over the sequences
  if (len(germ_seq) != len(seq)):
    print("Error the alignment sequence and the germline seq are with different sizes!!!")
    print(germ_seq)
    print(seq_name)
    print(seq)
    exit()
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
#IgPhyML output is masking the CDR3 region of the germline with gaps so we calculate the CDR3 consensus and insert it to Germline.
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
#Synonymous non Synonymous mutations counter for all sequences in clone lineage
#Function to count synonymous non synonymous mutations between germline sequence and SUT, sequence under test.
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
#Synonymous non Synonymous mutations counter for one sequence
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
#Synonymous non Synonymous mutations counter for one sequence with division to CDR and FR sections of the sequence
#Function to count synonymous non synonymous mutation ratio between germline sequence and farthest sequence and closest sequence to germline.
#Please note that in this function we calculate synonymous and not synonymous mutations for all amino acids.
def synonymous_mutation_ratio_per_seq_calc_cdr_fr_sections(germ_seq,sequence,cdr3_length):
  count_synon_cdr = 0
  count_not_synon_cdr = 0
  count_synon_fr = 0
  count_not_synon_fr = 0
  start_calc_flag = 0
  count_index = 1
  codon_str_seq = ""
  codon_str_germ = ""
  aa_str_seq = ""
  aa_str_germ = ""
  for index in range(0, len(sequence)):
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
            if (((index+1) > 79 and (index+1) <= 114) or ((index+1) > 166 and (index+1) <= 195) or ((index+1) > 310 and (index+1) <= (309+ cdr3_length))):  
                count_synon_cdr = count_synon_cdr + 1
            else:
                count_synon_fr = count_synon_fr + 1
        else:
            if (((index+1) > 79 and (index+1) <= 114) or ((index+1) > 166 and (index+1) <= 195) or ((index+1) > 310 and (index+1) <= (309+ cdr3_length))):
                count_not_synon_cdr = count_not_synon_cdr + 1
            else:
                count_not_synon_fr = count_not_synon_fr + 1
        codon_str_seq = ""
        codon_str_germ = ""
        count_index = 1
    else:
        count_index = count_index + 1
  if (count_synon_cdr == 0):
    sequence_omega_cdr = 0
  else:
    sequence_omega_cdr = (count_not_synon_cdr / count_synon_cdr)
  if (count_synon_fr == 0):
    sequence_omega_fr = 0
  else:
    sequence_omega_fr = (count_not_synon_fr / count_synon_fr)
  return sequence_omega_cdr, sequence_omega_fr

################################################################################
#Extraction of Sequence Alignment, Germline Alignment and CDR3 juntction from AIRR Table
#This extracts the sequence and the germline with the annotation from the AIRR DB file. this is required to calculate the CDR & FR Sections.
def airr_table_info_extraction(sequence_id, df_immune_db):
    id_original = sequence_id.replace("-",":")
    id_original = id_original.replace("M03592:154:000000000:","M03592:154:000000000-")
    airr_row_seq_id = df_immune_db.loc[df_immune_db['sequence_id'] == id_original]
    sequence_alignment = airr_row_seq_id.at[int(airr_row_seq_id.index.values),'sequence_alignment']
    germline_alignment = airr_row_seq_id.at[int(airr_row_seq_id.index.values),'germline_alignment']
    junction_cdr3 = airr_row_seq_id.at[int(airr_row_seq_id.index.values),'junction']
    return sequence_alignment, germline_alignment, len(junction_cdr3)
################################################################################
#Main Tree Analyzer Flow
#Tree Analyzer flow is the flow to extract all the tree features we are investigating per clone tree topology.
def tree_analyzer_flow(fasta,tree_topology_str,df_immune_db):
    # Load a tree and link it to an alignment from FASTA file.
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
    #Insert the annotated sequences to the tree:
    #Find the germline sequence and insert all other sequences to list:
    germline_sequence = ""
    germline_sequence_with_cdr3 = ""
    sequences_list = []
    for node in t.traverse(strategy="postorder"):
        temp = ""
        if node.is_leaf():
            if "GERM" in node.name:
                continue
            else:
                name_len = len(node.name)
                if (node.name[name_len-2:] != "_1"):
                    node.sequence, germline_sequence, temp = airr_table_info_extraction(node.name, df_immune_db)
                    sequences_list.append(node.sequence)
                    sequence_id_list.append(node.name)
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
################################################################################
    #Handle the germline CDR3 consensus:
    #germline_sequence_with_cdr3 = cdr3_consensus_calc(germline_sequence,sequences_list)
    germline_sequence_with_cdr3 = germline_sequence
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
            mutation_distance = mutation_dist_calc(germline_sequence_with_cdr3,node.sequence,node.name)
            if (mutation_distance >= max_distance_mutations):
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
                mutation_distance = mutation_dist_calc(nodeA.sequence,nodeB.sequence,nodeA.name)
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
    #Calculate non synonymous/synonymous mutations ratio for farthest and closest seq. to germline with division to CDR's and FR's sections
    max_sequence_alignment, max_germline_alignment, max_cdr3_length = airr_table_info_extraction(max_distance_seq_name, df_immune_db)
    omega_farthest_cdr, omega_farthest_fr = synonymous_mutation_ratio_per_seq_calc_cdr_fr_sections(max_germline_alignment,max_sequence_alignment,max_cdr3_length)
    min_sequence_alignment, min_germline_alignment, min_cdr3_length = airr_table_info_extraction(min_distance_seq_name, df_immune_db)
    omega_closest_cdr, omega_closest_fr = synonymous_mutation_ratio_per_seq_calc_cdr_fr_sections(min_germline_alignment,min_sequence_alignment,min_cdr3_length)
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
        "omega_ratio_farthest_sequence_cdr": omega_farthest_cdr,
        "omega_ratio_farthest_sequence_fr": omega_farthest_fr,
        "omega_ratio_closest_sequence_cdr": omega_closest_cdr,
        "omega_ratio_closest_sequence_fr": omega_closest_fr,
        "sequence_id_list": sequence_id_list
    }
    return results
