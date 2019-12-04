import jjm_analysis as jjm
import sys
import os


def normalize_in_directory(directory):


	for fname in os.listdir(directory):
	
		if fname.endswith('peaks.csv'): 

			array_name = fname[0:8] ;
		
			array_name = jjm.normalize_train_csv(fname) ; 

	return()
	
def normalize_in_directory_to_baseline(directory):

	#creates .csvs with normalized data
	#creates list of baseline files for further analysis
	
	baseline_file_list = []; 

	for fname in os.listdir(directory):
	
		if fname.endswith('peaks.csv'): 

			array_name = fname[0:8] ;
		
			array_name = jjm.normalize_train_csv_to_baseline(fname) ; 
			
			baseline_file = array_name; 

			baseline_file_list.append(baseline_file); 

	return(baseline_file_list)
			
			
			