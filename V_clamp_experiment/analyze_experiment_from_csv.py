import numpy as np
import csv
import pandas as pd
import os 
import sys
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts/V_clamp_experiment');
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/jupyterscripts');
import compile_by_abf_file

class analyzed_abf_file_pp40(object):
    ##class to store information for a stimfit .csv file
    
    def __init__(self, abf_file_name):
    #initialize file as name of able file
        self.abf_file_name = abf_file_name ; 
        
    def directory(self, root_dir):
    #returns the directory or directories containing the analyzed .csv file
        search_for = self.abf_file_name + '_stimfitanalysis.csv'
        found_dir = compile_by_abf_file.find_dir_iterative_for_extension(root_dir, search_for)[0] 
        return(found_dir)
		
class V_clamp_experiment(object):
    ## class to store data on voltage clamp/KA application experiments
    #object is a list of files with 1st file being the baseline
    
    def __init__(self, expt_name, expt_type, raw_files):
        #initialize object as a list of the raw data files
        self.raw_files = raw_files; 
        self.expt_name = expt_name;
        self.expt_type = expt_type; 
        
    def analyze_raw(self, raw_files, search_dir, baseline_indicies):
        ##here call     
        #cycle through raw_files, create analyzed file class to associate with directory
        print(search_dir)
        files_to_analyze = [];
        for fname in raw_files:
            stimfit_file = analyzed_abf_file_pp40(fname);
            files_to_analyze.append(str(stimfit_file.directory(search_dir))+'/'+str(fname)+'_stimfitanalysis.csv');
        print(files_to_analyze);
        os.chdir(stimfit_file.directory(search_dir)); 
        analyze_experiment_from_csv(files_to_analyze[0], baseline_indicies, str(self.expt_name), files_to_analyze[1:])
    
    def add_to_master(self, master_path, root_dir):
        ##updates master df with information by loading master, updating, and resaving
        #does not save
        path_to_master_df = '/Users/johnmarshall/Documents/Analysis/RecordingData/masterabfspreadsheet.xlsx';
        master_df_ = pd.read_excel(master_path);
        for abf_file in self.raw_files: 
            stimfit_file = analyzed_abf_file_pp40(abf_file); 
            d = {'associated experiment' : pd.Series([self.expt_name], index=['a']),
            'directory' : pd.Series([str(stimfit_file.directory(root_dir))], index=['a']),
            'experiment type': pd.Series([self.expt_type], index=['a']),
            'file': pd.Series([stimfit_file.abf_file_name], index=['a'])
            };
            to_add = pd.DataFrame(d);
            master_df_ = pd.concat([master_df_, to_add], axis = 0);
        return(master_df_)
    
    def update_master_file(self, path_to_master, path_to_save, root_dir):
        ##updates and saves master df 
        updated_master_df = self.add_to_master(path_to_master, root_dir); 
        updated_master_df.to_excel(path_to_save); 
        return(); 

def analyze_experiment_from_csv(baseline_csv_file, lines_from_end_for_baseline, expt_input, other_files):
	##analyze experiment function for reading csv files of analyzed data
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#get baseline from csv
	baseline_file_array = np.loadtxt(baseline_csv_file, delimiter=',', skiprows=1, usecols=[1,2,3,4,5,6,7])[-lines_from_end_for_baseline:];
	
	#mean values for baseline
	baseline_means = np.mean(baseline_file_array, axis = 0); 
	
	#divides arrays by mean values to create normalized values
	#does once for baseline array then variable times for subsequent depending on input
	baseline_normalized = np.divide(baseline_file_array, baseline_means); 
	
	experiment_timeseries_normalized = baseline_normalized; 
	
	for file in other_files:
		array = np.loadtxt(file, delimiter=',', skiprows=1, usecols=[1,2,3,4,5,6,7]);
		normalized = np.divide(array, baseline_means); 
		
	#concatenates with baseline 
		experiment_timeseries_normalized = np.vstack([experiment_timeseries_normalized, normalized]);  
	
	experiment_time_series_normalized = experiment_timeseries_normalized
	
	#make data frame with array
	column_names = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2']; 
	experiment_time_series_df = pd.DataFrame(experiment_time_series_normalized, columns=column_names);
	
	#make csv file with data
	expt_name = str(expt_input); 
	experiment_time_series_df.to_csv(expt_name + '_stimfitanalysis.csv', sep=',');
	
	
	return (experiment_time_series_normalized)