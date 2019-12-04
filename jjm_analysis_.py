import stf
from stf import get_base
import numpy as np
import csv


def jjm_resistance(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude):
	
	#time arguments in msec, amplitude argument in mV 
	
	stf.set_channel(0); 	
	stf.set_base_start(baseline_start, True) ; 
	stf.set_base_end(baseline_end, True) ; 
	stf.set_peak_start(cap_trans_start, True) ;
	stf.set_peak_end(cap_trans_end, True) ;
	stf.measure(); 
	
	baseline = float(stf.get_base()); 
	peak = float(stf.get_peak());
	
	real_peak = baseline - peak;
	
	amplitude = float(amplitude); 
	amplitude_V = amplitude/(10**(3)) ; 
	 
	real_peak_A = real_peak/(10**(12)) ;  
	
	Rs_Ohm = amplitude_V/abs(real_peak_A); 
	Rs = Rs_Ohm/(10**(6)) ; 
		
	return(real_peak, Rs)
	
	
def jjm_peak(baseline_start, baseline_end, p_start, p_end):

	#time arguments in msec, amplitude argument in mV 
	stf.set_channel(0); 	
	stf.set_base_start(baseline_start, True) ; 
	stf.set_base_end(baseline_end, True) ; 
	stf.set_peak_start(p_start, True) ;
	stf.set_peak_end(p_end, True) ;
	stf.measure();
	
	baseline = float(stf.get_base()); 
	peak = float(stf.get_peak());

	real_peak = abs(baseline - peak);

	return(real_peak)

	
def analyze_file(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end):

	"""inputs: (baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
	output: numpy array where 1st column is capacitance transient amplitude, 2nd is series resistance, 3rd is 1st EPSC, 4th is 2nd EPCSC
	also writes output to .csv file"""
	
	num_sweeps = stf.get_size_channel(); 
	print('there are')
	print(num_sweeps)
	print('sweeps in recording')
	print('analyzing sweeps')
	print(sweep_start)
	print('to')
	print(sweep_end)
	sweeps_to_analyze = sweep_end - sweep_start
	
	#create array for results
	data_array = np.zeros((sweeps_to_analyze+1, 4)) ; 

	y = 0 
	for x in range(sweep_start-1, sweep_end):
		#moves to next trace
		stf.set_trace(x);
		
		[cap_trans_amplitude, series_resistance] = jjm_resistance(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude);
		data_array[y][0] = cap_trans_amplitude;
		data_array[y][1] = series_resistance; 
		EPSC_1 = jjm_peak(baseline_start, baseline_end, EPSC1_s, EPSC1_e) ;
		data_array[y][2] = EPSC_1 ; 
		EPSC_2 = jjm_peak(baseline_start, baseline_end, EPSC2_s, EPSC2_e) ;
		data_array[y][3] = EPSC_2 ; 
		pp_40 = float(float(EPSC_2)/float(EPSC_1)); 	
		
		y += 1; 
	
	#print first few entries to check accuracy
	print(data_array[:3]); 
	
	#make csv file with data
	file_name = stf.get_filename(); 
	#expt = file_name[-12:].rstrip('.abf');
	np.savetxt(file_name + '_stimfitanalysis.csv', data_array , delimiter=',', newline='\n')
		
	return(data_array)
	
def import_test():
	print ('imported!');
	return()
	

def analyze_experiment(baseline_file_array, expt_input, *argv):
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#mean values for baseline
	baseline_means = np.mean(baseline_file_array, axis = 0); 
	
	#divides arrays by mean values to create normalized values
	#does once for baseline array then variable times for subsequent depending on input
	baseline_normalized = np.divide(baseline_file_array, baseline_means); 
	
	experiment_timeseries_normalized = baseline_normalized; 
	
	for array in argv:
		
		normalized = np.divide(array, baseline_means); 
		
	#concatenates with baseline 
		experiment_timeseries_normalized = np.vstack([experiment_timeseries_normalized, normalized]);  
	
	experiment_time_series_normalized = experiment_timeseries_normalized
	
	#create .csv file with concatenated data
	expt_name = str(expt_input); 
	np.savetxt(expt_name + '_stimfitanalysis.csv', experiment_time_series_normalized , delimiter=',', newline='\n'); 
	
	return (experiment_time_series_normalized)

def analyze_experiment_from_csv(baseline_csv_file, expt_input, *argv):
	##analyze experiment function for reading csv files of analyzed data
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#get baseline from csv
	baseline_file_array = np.loadtxt(baseline_csv_file, delimiter=',');
	
	#mean values for baseline
	baseline_means = np.mean(baseline_file_array, axis = 0); 
	
	#divides arrays by mean values to create normalized values
	#does once for baseline array then variable times for subsequent depending on input
	baseline_normalized = np.divide(baseline_file_array, baseline_means); 
	
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
	
def group_cells(sweep_groups_list, *argv):
	"first input is a list of tuples containing the indicies of sweeps to calculate values from in experiment:"
	"i.e. [(1,30),(31,45),(45,70)] 1:30 are baseline values, 31:45 are drug application, 45:70 wash"
	"next variable arguments are lists containing .csv files of normalized experiment data to group for analysis"
	"i.e. group 1: ['05062016_cell_1_stimfitanalysis.csv', '05062016_cell_2_stimfitanalysis.csv']"
	"output is a .csv file with group means for experimental values, should have first few columds with time series data," 
	"then following columns should list data from experiments" 
	
	data_by_group = []; 
	
	for experimental_group in argv:
		
		#here make an array of zeros for ultimate time series data
		#should be a 2d array with rows for the sum of the tuple input values
		#find length of experiment in sweeps
		length_of_experiment_periods = [];
		for x in range(len(sweep_groups_list)):
			length_of_experiment_periods.append(1+sweep_groups_list[x][1]-sweep_groups_list[x][0]);	
		print(length_of_experiment_periods);	
		#create array of zeros with 4 columns and total rows equalling length of experiment in sweeps
		experiment_length_in_sweeps = sum(length_of_experiment_periods);
		normalized_time_series_data = np.zeros((experiment_length_in_sweeps,4));

		#array for means of experiment periods
		experiment_means = np.zeros((experiment_length_in_sweeps,4));
		
		for file in experimental_group:
			
			#loads whole file from compiles time series .csv file appends to file with normalized data
			experiment_array = np.loadtxt(file, delimiter=',');
			experiment_sweeps_to_add = experiment_array[0:experiment_length_in_sweeps];
			normalized_time_series_data = np.hstack((normalized_time_series_data, experiment_sweeps_to_add));

			#adds normalized values to 1st 4 columns in array
			for column in range(4):
				for row in range(len(normalized_time_series_data)):
					normalized_time_series_data[row][column] = average_ith_value_in_array_row(normalized_time_series_data[row], 4, column+4, len(normalized_time_series_data[row]));
	
		data_by_group.append(normalized_time_series_data);
	
	return(normalized_time_series_data)

def calculate_means(sweep_groups_list, *argv):
	"similar to group_cells but calculates means of different periods in experiment"
	"first input is a list of tuples containing the indicies of sweeps to calculate values from in experiment:"
	"i.e. [(1,30),(31,45),(45,70)] 1:30 are baseline values, 31:45 are drug application, 45:70 wash"
	"next variable arguments are lists containing .csv files of normalized experiment data to group for analysis"
	means = [];
	for experimental_group in argv:
		group_means = [];
		for file in experimental_group:
			experiment_array = np.loadtxt(file, delimiter=',');
			
			#calculates means for values in experimental periods
			experiment_means = np.zeros((len(sweep_groups_list),4));
			for tuple in range(len(sweep_groups_list)):
				period_means = np.mean(experiment_array[(sweep_groups_list[tuple][0]-1):sweep_groups_list[tuple][1]], axis = 0);
				print(period_means);
				experiment_means[tuple] = period_means;	
			group_means.append(experiment_means);
		means.append(group_means);
	return(means)

			
def average_ith_value_in_array_row(array_row, i, start, stop):
	val = np.mean(array_row[start:stop:i]); 
	return(val)

def array_to_csv(array):
	np.savetxt('output.csv', array , delimiter=',', newline='\n');
	return()
		
def csv_to_array(csv_file):
	
	array = np.genfromtxt(csv_file, dtype=float, delimiter=',', names=True) ; 
	data_points = array['Input_0'] ; 
	
	time = array['Timems'] ; 
	sampling_interval = round(time[1]-time[0], 4);

	output_array = np.vstack((data_points, time)); 
	
	stf.new_window_matrix(output_array) ; 
	
	stf.set_sampling_interval(sampling_interval) ; 
	
	return(output_array)


		

def get_params(time):
	peak_1_s = stf.get_base_start(is_time=time);
	peak_1_e = stf.get_base_end(is_time=time);
	peak_2_s = stf.get_peak_start(is_time=time);
	peak_2_e = stf.get_peak_end(is_time=time);
	
	print(peak_1_s, peak_1_e, peak_2_s, peak_2_e); 
	params = [peak_1_s, peak_1_e, peak_2_s, peak_2_e]; 
	
	if time==True:
		print('values are time');
	else: 
		print('values are sweep indicies');
	
	return(params)
	
def set_params(params):	
	"""sets baseline and peak curors to input time values 
	(use timevalues not samples"""
	
	peak_1_s = params[0];
	peak_1_e = params[1];
	peak_2_s = params[2];
	peak_2_e = params[3];
	
	stf.set_base_start(peak_1_s, is_time=True);
	stf.set_base_end(peak_1_e, is_time=True);
	stf.set_peak_start(peak_2_s, is_time=True);
	stf.set_peak_end(peak_2_e, is_time=True);
	stf.measure();
	
	print(peak_1_s, peak_1_e, peak_2_s, peak_2_e); 
	
	return(params)
	
def increment_params(params, is_time, increment):
	"""increment cursors by input time points"""
		
	peak_1_s = params[0] + increment;
	peak_1_e = params[1] + increment;
	peak_2_s = params[2] + increment;
	peak_2_e = params[3] + increment;
	
	stf.set_base_start(peak_1_s, is_time=True);
	stf.set_base_end(peak_1_e, is_time=True);
	stf.set_peak_start(peak_2_s, is_time=True);
	stf.set_peak_end(peak_2_e, is_time=True);
	stf.measure();
	
	print(peak_1_s, peak_1_e, peak_2_s, peak_2_e); 
	
	return(params)

def increment_peak(params, is_time, increment):
	"""increment cursors by input time points"""
	
	peak_1_s = params[0];
	peak_1_e = params[1];		
	peak_2_s = params[2] + increment;
	peak_2_e = params[3] + increment;
	
	stf.set_base_start(peak_1_s, is_time=True);
	stf.set_base_end(peak_1_e, is_time=True);
	stf.set_peak_start(peak_2_s, is_time=True);
	stf.set_peak_end(peak_2_e, is_time=True);
	stf.measure();
	
	print(peak_1_s, peak_1_e, peak_2_s, peak_2_e); 
	
	params[2] = peak_2_s ; 
	params[3] = peak_2_e ;
	
	return(params)
	
def slice_peak_region(params, trace):
	
	"""use time for params, function converts to samples for cutting/displaying"""
	
	stf.select_trace(trace) ; 
	
	sampling_interval = stf.get_sampling_interval(); 
	
	peak_2_start_samples = (params[2] / sampling_interval) ;
	peak_2_end_samples = (params[3] / sampling_interval) ;
	
	peak_region = stf.get_trace()[peak_2_start_samples:peak_2_end_samples] ;  
	
	return(peak_region)

def scan_through_train(start_params, train_increment, num_stims, train_trace):
	
	"""scans through a tran of length "num_stims" in time increments of "train_increment", saves peak 
	amplitudes to an array peak_values (1st output) and sweep segments of peak regions for viewing are in peak_arrays"""
	
	stf.set_trace(train_trace) ; 
	
	baseline_s = start_params[0] ; 
	baseline_e = start_params[1] ; 
	params_ = start_params
	
	len_trace_in_samples = len(stf.get_trace(train_trace)); 
	peak_values = np.zeros(num_stims) ; 
	len_peak_region_in_samples = round((start_params[3] - start_params[2]) / stf.get_sampling_interval()) ; 
	peak_arrays = np.zeros((num_stims, (len_peak_region_in_samples))) ; 
	
	
	stim_count = 1; 
	while stim_count <= num_stims:
		 
		
		peak_start = params_[2] ;  
		peak_end = params_[3] ; 
		
		print(peak_start, peak_end) ;  
		peak = jjm_peak(baseline_s, baseline_e, peak_start, peak_end) ; 
		print(peak) ; 
		
		peak_values[stim_count-1] = peak ; 
	
		peak_region_slice = slice_peak_region(params_, train_trace) ; 
		
		peak_arrays[stim_count-1] = peak_region_slice ; 
		
		params_ = increment_peak(params_, True, train_increment) ; 
		
		stim_count += 1 ; 
	
	
	return(peak_values, peak_arrays)
	
def scan_through_train_expt(params_expt_input, train_increment, num_stims):
	
	
	
	len_peak_region_in_samples = round((params_expt_input[3] - params_expt_input[2]) / stf.get_sampling_interval())
	
	expt_peaks = np.zeros((stf.get_size_channel(), num_stims)) ; 
	expt_peak_arrays = np.zeros((stf.get_size_channel(), num_stims ,len_peak_region_in_samples)) ;
	
	trace = 0 ; 
	while trace < stf.get_size_channel():
		
		params_expt = params_expt_input ; 
		
		
		[expt_peaks[trace], expt_peak_arrays[trace]] = scan_through_train(params_expt, train_increment, num_stims, trace) ; 
		params_expt[2] = params_expt_input[2] - (train_increment*(num_stims)); 
		params_expt[3] = params_expt_input[3] - (train_increment*(num_stims));
		
		trace += 1; 
	
	loaded_file = stf.get_filename()[:-3] ; 
	np.savetxt(loaded_file + '_peaks.csv', expt_peaks , delimiter=',', newline='\n'); 
		
	return(expt_peaks, expt_peak_arrays) 
			
def batch_analysis(parameter_list, *argv):
	"inputs are list with analysis parameters(inputs for analyze file) and then *argv contains file names"
	
	batch_output = [];
	
	for file in argv:
		file_open(file);
		
		baseline_start = parameter_list[0]; 
		baseline_end = parameter_list[1]; 
		cap_trans_start = parameter_list[2];
		cap_trans_end = parameter_list[3];
		amplitude = parameter_list[4]; 
		EPSC1_s = parameter_list[5];
		EPSC1_e = parameter_list[6];
		EPSC2_s = parameter_list[7];
		EPSC2_e = parameter_list[8];
		sweep_start = parameter_list[9];
		sweep_end = parameter_list[10];
		
		file_array = analyze_file(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
		batch_output.append(file_array);
		
	return(batch_output)

def average_sweeps(*argv):
	
	sweeps = stf.get_trace(argv[0]);
	
	for sweep in argv[1:]:	
		sweep_ = stf.get_trace(sweep);
		sweeps = np.vstack((sweeps, sweep_)); 
	
	sweeps_mean = np.mean(sweeps, axis=0);
	
	stf.new_window(sweeps_mean); 
	
	return(sweeps_mean)
		
def remove_artifacts(art_1_start, art_1_end, art_2_start, art_2_end):
	sweep = stf.get_trace();
	artifacts_removed = np.hstack([sweep[:art_1_start],sweep[art_1_end:art_2_start],sweep[art_2_end:]]);
	return(artifacts_removed)






		
		