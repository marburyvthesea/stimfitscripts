import numpy as np

def analyze_experiment_from_csv(baseline_csv_file, expt_input, *argv):
	##analyze experiment function for reading csv files of analyzed data
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#get baseline from csv
	baseline_file_array = np.loadtxt(baseline_csv_file, delimiter=',');
	
	#mean values for baseline
	baseline_means = np.mean(baseline_file_array[-30:], axis = 0); 
	
	#divides arrays by mean values to create normalized values
	#does once for baseline array then variable times for subsequent depending on input
	baseline_normalized = np.divide(baseline_file_array[-30:], baseline_means); 
	
	experiment_timeseries_normalized = baseline_normalized; 
	
	for file in argv:
		array = np.loadtxt(file, delimiter=',');
		normalized = np.divide(array, baseline_means); 
		
	#concatenates with baseline 
		experiment_timeseries_normalized = np.vstack([experiment_timeseries_normalized, normalized]);  
	
	experiment_time_series_normalized = experiment_timeseries_normalized
	
	#create .csv file with concatenated data
	expt_name = str(expt_input); 
	np.savetxt(expt_name + '_stimfitanalysis.csv', experiment_time_series_normalized , delimiter=',', newline='\n'); 
	
	return (experiment_time_series_normalized)
	
def analyze_experiment_ppratio_from_csv(baseline_csv_file, expt_input, *argv):
	##analyze experiment pp ratio function for reading csv files of analyzed data
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#make list for ppratio time series
	#and means for separate experiment groups
	ppratio_timeseries = []; 
	ppratio_means = [];
	baseline_ppratio = [];
	
	#get baseline from csv
	baseline_file_array = np.loadtxt(baseline_csv_file, delimiter=',');
	
	#calculates ppratio by row and adds to lists	
	for row in baseline_file_array:
		ppratio_line = row[3]/row[2] ; 
		ppratio_timeseries.append(ppratio_line) ; 
		
			
	baseline_mean = np.mean(ppratio_timeseries, axis = 0); 
	ppratio_means.append(baseline_mean); 
	
	#goes through subsequent experiment files
	for file in argv:
		ppratio_experiment = []; 
		array = np.loadtxt(file, delimiter=',');
		for row in array:
			ppratio_line = row[3]/row[2] ; 
			ppratio_timeseries.append(ppratio_line);
			ppratio_experiment.append(ppratio_line);
			
		exp_mean = np.mean(ppratio_experiment, axis = 0);
		ppratio_means.append(exp_mean); 
	
	
	#create .csv file with data
	expt_name = str(expt_input); 
	np.savetxt(expt_name + '_ppratio_stimfitanalysis.csv', ppratio_timeseries , delimiter=',', newline='\n'); 
	
	return (ppratio_timeseries, ppratio_means)