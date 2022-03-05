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
def print_results_xsl(df):
	workbook = xlsxwriter.Workbook('Distribution_Analysis.xlsx')
	worksheet_spike_positive = workbook.add_worksheet('Spike+')
	worksheet_spike_negative = workbook.add_worksheet('Spike-')
	header_row = 1
	header_column = 0
	#Global Headers
	#spike negative
	for attribute_index in range(0,14):
		###############################################################################
		#Headers Rows
		worksheet_spike_negative.write(header_row,header_column,df.columns.values[attribute_index])
		worksheet_spike_negative.write(header_row+1,header_column,'Min')
		worksheet_spike_negative.write(header_row+2,header_column,'Max')
		worksheet_spike_negative.write(header_row+3,header_column,'Median')
		worksheet_spike_negative.write(header_row+4,header_column,'STD')
		#Headers Columns
		worksheet_spike_negative.write(header_row,header_column+1,'TP1')
		worksheet_spike_negative.write(header_row,header_column+2,'TP2')
		worksheet_spike_negative.write(header_row,header_column+3,'TP4')
		###############################################################################
		for meas_index in range(0,4):
			for tp_index in range(0,3):
				worksheet_spike_negative.write(header_row+meas_index+1,header_column+tp_index+1,results_container[0,meas_index,tp_index,attribute_index])
		#Update the header row for next attribute
		header_row += 6
	#spike positive
	header_row = 0
	header_column = 0
	for attribute_index in range(0,14):
		###############################################################################
		#Headers Rows
		worksheet_spike_positive.write(header_row,header_column,df.columns.values[attribute_index])
		worksheet_spike_positive.write(header_row+1,header_column,'Min')
		worksheet_spike_positive.write(header_row+2,header_column,'Max')
		worksheet_spike_positive.write(header_row+3,header_column,'Median')
		worksheet_spike_positive.write(header_row+4,header_column,'STD')
		#Headers Columns
		worksheet_spike_positive.write(header_row,header_column+1,'TP1')
		worksheet_spike_positive.write(header_row,header_column+2,'TP2')
		worksheet_spike_positive.write(header_row,header_column+3,'TP4')
		###############################################################################
		for meas_index in range(0,4):
			for tp_index in range(0,3):
				worksheet_spike_positive.write(header_row+meas_index+1,header_column+tp_index+1,results_container[1,meas_index,tp_index,attribute_index])
		#Update the header row for next attribute
		header_row += 6
	workbook.close()
################################################################################
#Read Tree Analyzer Merged Results
tree_analyzer_df = pd.read_csv("subject022_test.csv")
#Divide to Spike+ and Spike- Results
spike_positive_df = tree_analyzer_df.loc[tree_analyzer_df['spike'] == '+']
spike_negative_df = tree_analyzer_df.loc[tree_analyzer_df['spike'] == '-']
analyze_distribution(spike_negative_df,0)
analyze_distribution(spike_positive_df,1)
np.set_printoptions(suppress=True)
#print(results_container)
print_results_xsl(tree_analyzer_df)

