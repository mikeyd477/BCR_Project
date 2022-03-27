"""

author: Michael Dubovsky

moodule to analyzer statistical attributes of the tree analyzer distribution

"""
################################################################################
#Imports
import pandas as pd
import time
import numpy as np
import xlsxwriter
import argparse
import os
################################################################################
#Start Time
start_time = time.time()
#Array to hold results spike positve index is 1, spike negative index is 0
results_container = np.zeros((2,4,3,14)) 
################################################################################
#Function to analyze statistical qunatities of the results attributes
#Calculates Min,Max,Median,Standart Deviation.
def analyze_distribution(df,spike_indicator):
	tp_indicator = 0
	for i in range(0,3):
		if (i == 2):
			tp_indicator = 4
		else:
			tp_indicator = i+1
		tp_df = df.loc[df['TP'] == tp_indicator]
		tp_df_min = tp_df.min()
		tp_df_max = tp_df.max()
		tp_df_median = tp_df.median()
		tp_df_std = tp_df.std()
		for j in range(0,14):
			results_container[spike_indicator,0,i,j] = tp_df_min[j]
			results_container[spike_indicator,1,i,j] = tp_df_max[j]
			results_container[spike_indicator,2,i,j] = tp_df_median[j]
			results_container[spike_indicator,3,i,j] = tp_df_std[j]
################################################################################
#Function to Write Results to EXCEL XLS FILE
#
def print_results_xsl(df,output):
	workbook = xlsxwriter.Workbook(output)
	worksheet = workbook.add_worksheet('Statistical_Summary')
	header_row = 1
	header_column = 0
	#Global Headers
	worksheet.write(0,2,'Spike-')
	worksheet.write(0,8,'Spike+')
	#spike negative
	for attribute_index in range(0,14):
		###############################################################################
		#Headers Rows
		worksheet.write(header_row,header_column,df.columns.values[attribute_index])
		worksheet.write(header_row+1,header_column,'Min')
		worksheet.write(header_row+2,header_column,'Max')
		worksheet.write(header_row+3,header_column,'Median')
		worksheet.write(header_row+4,header_column,'STD')
		#Headers Columns
		worksheet.write(header_row,header_column+1,'TP1')
		worksheet.write(header_row,header_column+2,'TP2')
		worksheet.write(header_row,header_column+3,'TP4')
		###############################################################################
		for meas_index in range(0,4):
			for tp_index in range(0,3):
				worksheet.write(header_row+meas_index+1,header_column+tp_index+1,results_container[0,meas_index,tp_index,attribute_index])
		#Update the header row for next attribute
		header_row += 6
	#spike positive
	header_row = 1
	header_column = 6
	for attribute_index in range(0,14):
		###############################################################################
		#Headers Rows
		worksheet.write(header_row,header_column,df.columns.values[attribute_index])
		worksheet.write(header_row+1,header_column,'Min')
		worksheet.write(header_row+2,header_column,'Max')
		worksheet.write(header_row+3,header_column,'Median')
		worksheet.write(header_row+4,header_column,'STD')
		#Headers Columns
		worksheet.write(header_row,header_column+1,'TP1')
		worksheet.write(header_row,header_column+2,'TP2')
		worksheet.write(header_row,header_column+3,'TP4')
		###############################################################################
		for meas_index in range(0,4):
			for tp_index in range(0,3):
				worksheet.write(header_row+meas_index+1,header_column+tp_index+1,results_container[1,meas_index,tp_index,attribute_index])
		#Update the header row for next attribute
		header_row += 6
	workbook.close()
################################################################################
#Parse Command Line Argument Section
#Create object for parsing command-line options
parser = argparse.ArgumentParser(description="Command line for Distribution Analyzer")
#Add arguments to command line parser
parser.add_argument("-d", "--data", type=str, help="Path to data of distribution output from tree parser")
parser.add_argument("-o", "--output", type=str, help="Output File Name")
#Parse the command line arguments to an object
args = parser.parse_args()
#Safety if no parameter have been given
if not args.data:
	print("No Data File has been given. Please insert DATA file in command line --data")
	print("For help type --help")
	exit()
if not args.output:
	print("No Output File Name has been given. Please insert Output File Name in command line --output")
	print("For help type --help")
	exit()
################################################################################
#Read Tree Analyzer Merged Results
tree_analyzer_df = pd.read_csv(args.data)
#Divide to Spike+ and Spike- Results
spike_positive_df = tree_analyzer_df.loc[tree_analyzer_df['spike'] == '+']
spike_negative_df = tree_analyzer_df.loc[tree_analyzer_df['spike'] == '-']
analyze_distribution(spike_negative_df,0)
analyze_distribution(spike_positive_df,1)
np.set_printoptions(suppress=True)
#print(results_container)
print_results_xsl(tree_analyzer_df,args.output)

