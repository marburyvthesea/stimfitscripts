import numpy as np
import csv
import os
import sys
import math

add_to_path = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/wx-3.0-osx_cocoa', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks/stimfit', '.', '', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Applications/MacPorts/stimfit.app/Contents/Frameworks/stimfit', '/Applications/MacPorts/stimfit.app/Contents/Frameworks', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Library/Python/2.7/site-packages']
more_add_to_path = ['/anaconda3/envs/py27/lib/python2.7/site-packages','', '/Library/Python/2.7/site-packages/openpyxl-2.3.0-py2.7.egg', '/Library/Python/2.7/site-packages/et_xmlfile-1.0.1-py2.7.egg', '/Library/Python/2.7/site-packages/jdcal-1.2-py2.7.egg', '/Library/Python/2.7/site-packages/XlsxWriter-0.7.7-py2.7.egg', '/Library/Python/2.7/site-packages', '/Library/Python/2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC', '/Library/Python/2.7/site-packages'];
for path in add_to_path:
	sys.path.append(path);  
for path in more_add_to_path:
	sys.path.append(path);
import pandas as pd	
from pandas import ExcelWriter
import stf
from stf import get_base
import convert_prairie_tseries_to_abf
import prairie_tseries_toabf_stimfit as pv



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
	
	try:
		Rs_Ohm = amplitude_V/abs(real_peak_A); 
		Rs = Rs_Ohm/(10**(6)) ; 
	except ZeroDivisionError:
		print('Cap trans peak not found'); 
		Rs = 'NaN' ; 
		pass 
				
	return(real_peak, Rs)

def jjm_resistance_(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude):
	
	#time arguments in msec, amplitude argument in mV 
	#array for output
	np.array((stf.get_size_channel(), 2), dtype=[('cap trans peak(pA)', '<f4'), ('Series Resistance', '<f4')]);  
	
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
	
def jjm_mean(baseline_start, baseline_end, m_start, m_end):

	#time arguments in msec, amplitude argument in mV 
	stf.set_channel(0); 	
	stf.set_base_start(baseline_start, True) ; 
	stf.set_base_end(baseline_end, True) ; 
	
	stf.set_peak_start(m_start, True) ; 
	stf.set_peak_end(m_end, True) ; 
	
	
	baseline = float(stf.get_base()); 
	
	sample_point_start = int(round(m_start/stf.get_sampling_interval())); 
	sample_point_end = int(round(m_end/stf.get_sampling_interval())); 
	
	mean_region = stf.get_trace()[sample_point_start:sample_point_end];
	mean = np.mean(mean_region); 
	
	real_mean = abs(baseline - mean);

	
	return(real_mean)


	
def analyze_file(baseline_start, baseline_end, baseline_2_start, baseline_2_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end):

	"""inputs: (baseline_start, baseline_end, baseline_2_start, baseline_2_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
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
	data_array = np.zeros((sweeps_to_analyze, 7)) ; 

	y = 0 
	for x in range(sweep_start, sweep_end):
		#moves to next trace
		stf.set_trace(x);
		
		#jjm_resistance to get series resistance values 
		[cap_trans_amplitude, series_resistance] = jjm_resistance(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude);
		data_array[y][0] = cap_trans_amplitude;
		data_array[y][1] = series_resistance; 
		
		#input resistance calculation
		[SS, Ri] = input_resistance(cap_trans_end, amplitude);
		data_array[y][2] = SS; 
		data_array[y][3] = Ri; 
		
		#jjm_capacitance to get whole cell capacitance estimate
		whole_cell_capacitance = capacitance_reimann([baseline_start, baseline_end, cap_trans_start, cap_trans_end], amplitude);
		data_array[y][4] = whole_cell_capacitance; 
		
		#EPSC values
		EPSC_1 = jjm_peak(baseline_start, baseline_end, EPSC1_s, EPSC1_e) ;
		data_array[y][5] = EPSC_1 ; 
		EPSC_2 = jjm_peak(baseline_2_start, baseline_2_end, EPSC2_s, EPSC2_e) ;
		data_array[y][6] = EPSC_2 ; 
		try:
			pp_40 = float(float(EPSC_2)/float(EPSC_1)); 	
		except ZeroDivisionError:
			print('could not find EPSC peaks'); 
			pp_40 = 'NaN' ; 
		
		y += 1; 
	
	#print first few entries to check accuracy
	print(data_array[:3]); 
	
	#make data frame with array
	column_names = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2']; 
	data_array_df = pd.DataFrame(data_array, columns=column_names);
	
	#make csv file with data
	file_name = stf.get_filename(); 
	data_array_df.to_csv(file_name + '_stimfitanalysis.csv', sep=',')
		
	return(data_array_df)
	
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

def analyze_experiment_from_csv(path_to_master_, experiment_type_, expt_input, baseline_csv_file, *argv):
	##analyze experiment function for reading csv files of analyzed data
	##want this function to take numpy array for baseline values, average all columns to create average baseline values
	#then normalize all the experiment numpy arrays to the baseline values to create one numpy time series for the experiment
	
	#get baseline from csv
	baseline_file_array = np.loadtxt(baseline_csv_file, delimiter=',', skiprows=1, usecols=[1,2,3,4,5,6,7]);
	
	#mean values for baseline
	baseline_means = np.mean(baseline_file_array[-31:], axis = 0); 
	
	#divides arrays by mean values to create normalized values
	#does once for baseline array then variable times for subsequent depending on input
	baseline_normalized = np.divide(baseline_file_array[-31:], baseline_means); 
	
	experiment_timeseries_normalized = baseline_normalized; 
	
	for f in argv:
		array = np.loadtxt(f, delimiter=',', skiprows=1, usecols=[1,2,3,4,5,6,7]);
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
	
	expt_files = argv
	analyze_experiment_update_master(path_to_master_, experiment_type_, baseline_csv_file, expt_input, expt_files)
	
	return (experiment_time_series_normalized)

def analyze_experiment_update_master(path_to_master, experiment_type, baseline_csv_file_, expt_input_, *argv):
	"""Normalizes csv files to an input baseline csv file, updates master abf sheet(at location indicated) to associate raw abf files with
	an experiment file"""
	
	abf_files = list(argv[0]);
	master_abf_df = pd.read_excel(path_to_master, delimiter=',');
	abf_files.append(baseline_csv_file_);
	print(abf_files)
	for fname in abf_files:
		print(fname);
		files_list_df = master_abf_df.file.as_matrix();
		#adds file to df if not in already
		if fname not in files_list_df:
			to_add = {u'file': str(fname[:12])};
			print('adding', to_add['file']); 
			master_abf_df = master_abf_df.append(to_add, ignore_index=True);
		#updates sheet with experiment information
		master_abf_df.loc[master_abf_df.file==fname[:12], ['associated experiment', 'experiment type']]= expt_input_, experiment_type ; 
	

	master_abf_df.to_excel(path_to_master);
	
	
	#analyze_experiment_from_csv(baseline_csv_file_, expt_input_, abf_files); 	
		
	return()

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

def array_to_named_csv(array, filename):
	np.savetxt(filename+'.csv', array , delimiter=',', newline='\n');
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
	
	baseline_s = params[0] + increment
	baseline_e = params[1] + increment	
	peak_2_s = params[2] + increment
	peak_2_e = params[3] + increment
	
	stf.set_base_start(baseline_s, is_time=True)
	stf.set_base_end(baseline_e, is_time=True)
	stf.set_peak_start(peak_2_s, is_time=True)
	stf.set_peak_end(peak_2_e, is_time=True)
	stf.measure()
	
	print(baseline_s, baseline_e, peak_2_s, peak_2_e)
	
	params[0] = baseline_s 
	params[1] = baseline_e
	params[2] = peak_2_s 
	params[3] = peak_2_e 
	
	return(params)
	
def slice_peak_region(params, trace):
	
	"""use time for params, function converts to samples for cutting/displaying"""
	
	stf.select_trace(trace) ; 
	
	sampling_interval = stf.get_sampling_interval(); 
	
	peak_2_start_samples = int(params[2] / sampling_interval) ;
	peak_2_end_samples = int(params[3] / sampling_interval) ;
	
	peak_region = stf.get_trace()[peak_2_start_samples:peak_2_end_samples] ;  
	
	return(peak_region)

def scan_through_train(start_params, train_increment, num_stims, train_trace):
	
	"""scans through a tran of length "num_stims" in time increments of "train_increment", saves peak 
	amplitudes to an array peak_values (1st output) and sweep segments of peak regions for viewing are in peak_arrays"""
	
	stf.set_trace(train_trace) 
	
	params_ = start_params
	
	len_trace_in_samples = int(len(stf.get_trace(train_trace)))
	peak_values = np.zeros(num_stims) 
	len_peak_region_in_samples = round((start_params[3] - start_params[2]) / stf.get_sampling_interval()) 
	peak_arrays = np.zeros((int(num_stims), int(len_peak_region_in_samples))) 
	
	stim_count = 1
	while stim_count <= num_stims:
		 
		print(params_[0], params_[1], params_[2], params_[3]) 
		peak = jjm_peak(params_[0], params_[1], params_[2], params_[3]) 
		print(peak)  
		
		peak_values[int(stim_count-1)] = peak 
	
		peak_region_slice = slice_peak_region(params_, train_trace) 
		
		if len(peak_arrays[int(stim_count-1)]) == len(peak_region_slice):
			peak_arrays[int(stim_count-1)] = peak_region_slice
		
		elif len(peak_arrays[int(stim_count-1)]) > len(peak_region_slice):
			to_pad = len(peak_arrays[int(stim_count-1)]) - len(peak_region_slice)
			peak_region_slice = np.pad(peak_region_slice, (0, to_pad), 'constant', constant_values=(0, to_pad))
			peak_arrays[int(stim_count-1)] = peak_region_slice
			
		else: 
			peak_arrays[int(stim_count-1)] = peak_region_slice[0:len(peak_arrays[int(stim_count-1)])]
			
		# increments peak and baseline values by train stim frequency
		params_ = increment_params(params_, True, train_increment) 
		
		stim_count += 1 
	
	return(peak_values, peak_arrays)
	
def scan_through_train_expt(params_expt_input, train_increment, num_stims):
	
	len_peak_region_in_samples = round((params_expt_input[3] - params_expt_input[2]) / stf.get_sampling_interval())
	
	expt_peaks = np.zeros((int(stf.get_size_channel()+1), int(num_stims)), dtype='float32') ; 
	expt_peak_arrays = np.zeros((int(stf.get_size_channel()), int(num_stims) , int(len_peak_region_in_samples))) ;
	
	trace = 0 ; 
	while trace < stf.get_size_channel():
		
		print('params should reset to')
		print(params_expt_input)
		params_expt = params_expt_input ; 
		
		[expt_peaks[trace], expt_peak_arrays[trace]] = scan_through_train(params_expt, train_increment, num_stims, trace) ; 
		params_expt[0] = params_expt_input[0] - (train_increment*(num_stims)); 
		params_expt[1] = params_expt_input[1] - (train_increment*(num_stims));
		params_expt[2] = params_expt_input[2] - (train_increment*(num_stims)); 
		params_expt[3] = params_expt_input[3] - (train_increment*(num_stims));
		
		trace += 1; 
	
	#expt_peaks[stf.get_size_channel()+1, 0] = 'file params'
	#expt_peaks[stf.get_size_channel()+1, 1] = 'file params'
	
	loaded_file = stf.get_filename()[:-3] ; 
	np.savetxt(loaded_file + '_all_peaks.csv', expt_peaks , delimiter=',', newline='\n'); 
		
	return(expt_peaks, expt_peak_arrays) 
	
def scan_pp40(pp40_params):
	
	pp40_output_array_peaks = scan_through_train_expt(pp40_params, 40, 2)[0] ; 
			
	return(pp40_output_array_peaks)
	
def pp40_expt(pp40_params, captrans_params):
	
	##make this calc captrans for each sweepx
	cap_trans_output_peaks = jjm_resistance(captrans_params[0], captrans_params[1], captrans_params[2], captrans_params[3], captrans_params[4]) ; 	
	pp40_output_array_peaks = scan_through_train_expt(pp40_params, 40, 2)[0] ; 
	pp40_expt_output = np.vstack((cap_trans_output_peaks, pp40_output_array_peaks));   
	
	loaded_file = stf.get_filename()[:-3] ; 
	np.savetxt(loaded_file + '_pp40_peaks.csv', pp40_expt_output  , delimiter=',', header="EPSC Amplitude (pA)", newline='\n');
	
	return(pp40_expt_output, cap_trans_output_peaks )
	
def batch_analysis(parameter_list, file_list):
	"""inputs are list with analysis parameters(inputs for analyze file) and then *argv contains file names"
	analysis parameters: (baseline_2_start, baseline_2_end, baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, 
	EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweeps_from_end_to_analyze or 'all'
	"""
	
	batch_output_names = []
	batch_output_dfs = []
	batch_output_dict = {}
	
	for file in file_list:
		print('file is:')
		print(file)
		if file.endswith('.abf'):
			stf.file_open(file);
		else:
			if file.startswith('TSeries'):
				sweeps_compiled_from_pv_tseries = pv.import_t_series_episodic(file);
				print('loaded_directory')
				pv.plot_episodic_array(sweeps_compiled_from_pv_tseries); 
				
		baseline_start = int(parameter_list[0]); 
		baseline_end = int(parameter_list[1]); 
		baseline_2_start = int(parameter_list[2]); 
		baseline_2_end = int(parameter_list[3]); 
		cap_trans_start = int(parameter_list[4]);
		cap_trans_end = int(parameter_list[5]);
		amplitude = int(parameter_list[6]); 
		EPSC1_s = int(parameter_list[7]);
		EPSC1_e = int(parameter_list[8]);
		EPSC2_s = int(parameter_list[9]);
		EPSC2_e = int(parameter_list[10]);
		# either analyze a set number of sweeps from end or if 'all' then analyze all sweeps
		if parameter_list[11] == 'all':
			sweep_start = 0 ;
			sweep_end = stf.get_size_channel() ;
		else:
			sweeps_from_end_to_analyze = int(parameter_list[11]); 
			sweep_start = stf.get_size_channel()-sweeps_from_end_to_analyze;
			sweep_end = stf.get_size_channel();
		
		data_array = analyze_file(baseline_start, baseline_end, baseline_2_start, baseline_2_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
		#calculate pp40 column
		
		data_array_df = pd.concat([data_array, pd.DataFrame({'pp40':data_array['EPSC2']/data_array['EPSC1']})], axis =1 )

	
		#make csv file with data
		data_array_df.to_csv(file + '_stimfitanalysis.csv', sep=',');
		batch_output_dict[str(file)] = data_array_df
		batch_output_dfs.append(data_array_df)
		batch_output_names.append(file)
		
	batch_output = pd.concat(batch_output_dict)
		
	return(batch_output)


def pp40_expt_analysis(parameter_list, *argv):
	"inputs are list with analysis parameters(inputs for analyze file) and then *argv contains file names"
	
	batch_output = {};
	
	file_count = 0 ; 
	while file_count < len(argv):
		if argv[file_count].endswith('.abf'):
			stf.file_open(argv[file_count]);
		else:
			if file.startswith('TSeries'):
				sweeps_compiled_from_pv_tseries = pv.import_t_series_episodic(argv[file_count]);
				pv.plot_episodic_array(sweeps_compiled_from_pv_tseries); 
		
		baseline_start = parameter_list[0]; 
		baseline_end = parameter_list[1]; 
		cap_trans_start = parameter_list[2];
		cap_trans_end = parameter_list[3];
		amplitude = parameter_list[4]; 
		EPSC1_s = parameter_list[5];
		EPSC1_e = parameter_list[6];
		EPSC2_s = parameter_list[7];
		EPSC2_e = parameter_list[8];
		sweeps_from_end_to_analyze = int(parameter_list[9][file_count]); 
		sweep_start = stf.get_size_channel()-sweeps_from_end_to_analyze;
		sweep_end = stf.get_size_channel();
		
		file_array = analyze_file(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
		array_to_named_csv(file_array, argv[file_count]); 
		batch_output[argv[file_count]] = file_array ;
	
		file_count += 1; 	
	return(batch_output)



def NMDA_AMPA_analysis(parameter_list, *argv):
	"""inputs are parameters, then files to analyze
	inputs: (baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, NMDA_threshold_search_start, NMDA_threshold_search_end)
	NMDA measurements are the regions in which to search for the NMDA stim artifact
	baseline is then roughly 10ms before artifact, NMDA amplitude is mean value over 2ms, 60ms from artifact"""
	NMDA_AMPA_batch_output = {};
	
	for file in argv:
		##checks if abf or prairieview tseries and handles accordingly
		if file.endswith('.abf'):
			stf.file_open(file);
		else:
			if file.startswith('TSeries'):
				sweeps_compiled_from_pv_tseries = pv.import_t_series_episodic(file);
				pv.plot_episodic_array(sweeps_compiled_from_pv_tseries); 
		
		baseline_start = parameter_list[0]; 
		baseline_end = parameter_list[1]; 
		cap_trans_start = parameter_list[2];
		cap_trans_end = parameter_list[3];
		amplitude = parameter_list[4]; 
		EPSC1_s = parameter_list[5];
		EPSC1_e = parameter_list[6];
		EPSC2_s = EPSC1_s;
		EPSC2_e = EPSC1_e;
		##will analyze all sweeps in file
		sweep_start = 0;
		sweep_end = stf.get_size_channel();
		
		##analyzes AMPA values, for now ignore 2nd EPSC call
		file_array = analyze_file(baseline_start, baseline_end, baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end);
		
		##scans through sweeps and finds NMDA peaks
		stf.set_trace(0); 
		NMDA_peaks = []
		for sweep in range(sweep_end):
			print(sweep)
			stf.set_trace(sweep); 
			NMDA_threshold_search_start = parameter_list[7];
			NMDA_threshold_search_end = parameter_list[8];
		
			threshold_time = jjm_count(NMDA_threshold_search_start, NMDA_threshold_search_end-NMDA_threshold_search_start, 400, up=True, trace=sweep, mark=True)[2][0]
		
			stf.set_trace(sweep); 
			
			NMDA_baseline_start = threshold_time-13;
			
			NMDA_baseline_end = threshold_time-3;
			
			NMDA_EPSC1_s = threshold_time+60;	
			
			NMDA_EPSC1_e = threshold_time+62;
				
			##calls jjm_mean
			stf.set_trace(sweep);
			set_params([NMDA_baseline_start, NMDA_baseline_end, NMDA_EPSC1_s, NMDA_EPSC1_e]) 
			stf.measure()
			NMDA_peak = jjm_mean(NMDA_baseline_start, NMDA_baseline_end, NMDA_EPSC1_s, NMDA_EPSC1_e); 	
					
			##paste in NMDA values into AMPA value array column 3
			NMDA_peaks.append(NMDA_peak) 
		
		files_compiled = pd.concat([file_array, pd.DataFrame({'NMDA_peaks': NMDA_peaks})], axis=1)
		
		files_compiled.to_excel(file + 'AMPA_NMDA_stimfitanalysis.xlsx')
		
	return(NMDA_peaks, files_compiled)
	
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

def remove_artifacts_untrim(art_1_start, art_1_end, art_2_start, art_2_end):
	sweep = stf.get_trace();
	artifacts_removed = np.hstack([sweep[:art_1_start],sweep[art_1_end:art_2_start],sweep[art_2_end:]]);
	return(artifacts_removed)


def normalize_train_csv(csv_file): 
	try:
		array = np.loadtxt(csv_file, dtype=float, delimiter=',') ;
		_1st_EPSC = array[0][0]; 
		norm_array = array/_1st_EPSC
		expt_name = csv_file.rstrip('.csv')
		np.savetxt(expt_name + '_normalized.csv', np.transpose(norm_array) , delimiter=',', newline='\n'); 
	except IndexError:
		print('IndexError')
		print(csv_file); 
		pass
	except ValueError:
		print('ValueError')
		print(csv_file)
		pass
	return() 

def normalize_train_csv_to_baseline_input(csv_file, baseline_file): 
	array = np.loadtxt(csv_file, dtype=float, delimiter=',') ;
	baseline_mean_1st_EPSC = np.mean(np.transpose(np.loadtxt(baseline_file, dtype=float, delimiter=','))[2]); 
	norm_array = array/baseline_mean_1st_EPSC 
	expt_name = csv_file.rstrip('.csv')
	np.savetxt(expt_name + '_normalized_baseline.csv', np.transpose(norm_array) , delimiter=',', newline='\n'); 
	return()
	
def normalize_train_csv_to_baseline(csv_file): 
	
	## returns baseline file for further analysis
	
	try:
		array = np.loadtxt(csv_file, dtype=float, delimiter=',') ;
		file_index = int(csv_file[6:8]); 
		baseline_index = file_index-1; 
		baseline_index_str = str(baseline_index); 
		
		if len(baseline_index_str) < 2:
			baseline_index_str = '0' + baseline_index_str; 
		
		baseline_file = csv_file[:6] + baseline_index_str + '.abf_stimfitanalysis.csv' ; 
		
		print('loading', baseline_file)
		## first tries to find a file by just adding '.abf_stimfitanalysis.csv' to the previous index
		try:
			baseline_mean_1st_EPSC = np.mean(np.transpose(np.loadtxt(baseline_file, dtype=float, delimiter=',', skiprows=1))[6][-15:])
			print('mean is', baseline_mean_1st_EPSC)
			norm_array = array/baseline_mean_1st_EPSC 
			expt_name = csv_file.rstrip('.csv')
			np.savetxt(expt_name + '_normalized_baseline.csv', np.transpose(norm_array) , delimiter=',', newline='\n'); 
		except IOError:
			
			## checks if baseline file contains 'filtered' string	
			try: 
				baseline_file_filtered = csv_file[:6] + baseline_index_str + '_filter1k.abf_stimfitanalysis.csv' ;	
				baseline_mean_1st_EPSC = np.mean(np.transpose(np.loadtxt(baseline_file_filtered, dtype=float, delimiter=','))[2])
				baseline_file = baseline_file_filtered ; 
				norm_array = array/baseline_mean_1st_EPSC 
				expt_name = csv_file.rstrip('.csv')
				np.savetxt(expt_name + '_normalized_baseline.csv', np.transpose(norm_array) , delimiter=',', newline='\n'); 
	
			except IOError:
				print('IOError: could not find')
				print(baseline_file)
				print('or')
				print(baseline_file_filtered)
				print('baseline for')
				print(csv_file)
				baseline_file = 'none';
				pass
				
			except IndexError:
				print('IndexError')
				print(csv_file); 
				baseline_file = 'none';
				pass
			
			except ValueError:
				print('ValueError')
				print(csv_file)
				baseline_file = 'none';
				pass
	
	except IndexError:
		print('IndexError')
		print(csv_file); 
		baseline_file = 'none'; 	
		pass
	
	except ValueError:
		print('ValueError')
		print(csv_file)
		baseline_file = 'none';
		pass

	return(baseline_file) 

def compile_1st_row(file_list):

	compiled_array = np.zeros((len(file_list), 5)); 
	file_count = 0 ; 
	while file_count < len(file_list):
		file = file_list[file_count]; 
		array = np.loadtxt(file, delimiter=',');
		compiled_array[file_count] = array[0] ;
		file_count += 1 ; 	
		
	return(compiled_array)


def normalize_in_directory(directory):
	for fname in os.listdir(directory):
		if fname.endswith('peaks.csv'): 
			try:	
				normalize_train_csv(fname) ; 
			except ValueError:
				pass
	return()
	
def jjm_count(start, delta, threshold=0, up=True, trace=None, mark=True):
	""" Counts the number of events (e.g action potentials (AP)) in the current trace.
	Arguments:
	start -- starting time (in ms) to look for events.
	delta -- time interval (in ms) to look for events.
	threshold -- (optional) detection threshold (default = 0).
	up -- (optional) True (default) will look for upward events, False downwards.
	trace -- (optional) zero-based index of the trace in the current channel,
	if None, the current trace is selected.
	mark -- (optional) if True (default), set a mark at the point of threshold crossing
	Returns:
	An integer with the number of events.
	Examples:
	count_events(500,1000) returns the number of events found between t=500 ms and t=1500 ms
	above 0 in the current trace and shows a stf marker.
	count_events(500,1000,0,False,-10,i) returns the number of events found below -10 in the trace i and shows the corresponding stf markers. """
	# sets the current trace or the one given in trace.
	if trace is None:
		sweep = stf.get_trace_index()
	else:
		if type(trace) !=int:
			print "trace argument admits only integers"
			return False
		sweep = trace
	# set the trace described in sweep
	stf.set_trace(sweep)
	# transform time into sampling points
	dt = stf.get_sampling_interval()
	pstart = int( round(start/dt) )
	pdelta = int( round(delta/dt) )
	# select the section of interest within the trace
	selection = stf.get_trace()[pstart:(pstart+pdelta)]
	# algorithm to detect events
	EventCounter,i = 0,0 # set counter and index to zero
	# list of sample points
	sample_points = []
	# choose comparator according to direction:
	if up:
		comp = lambda a, b: a > b
	else:
		comp = lambda a, b: a < b
	# run the loop
	while i<len(selection):
		if comp(selection[i],threshold):
			EventCounter +=1
			if mark:
				sample_point = pstart+i; 
				sample_points.append(sample_point); 
				stf.set_marker(pstart+i, selection[i])
			while i<len(selection) and comp(selection[i],threshold):
				i+=1 # skip values if index in bounds AND until the value is below/above threshold again
		else:
			i+=1
	
	time_points = [];
	for point in sample_points:
		time_points.append(sample_point*dt); 
	return EventCounter, sample_points, time_points

def jjm_fit(params, fit_end_time, function_fit):
	##TO DO: find only positive going peak
	"""returns function fit to fit cursors
		Inputs: params, list of [baseline_start, baseline_end, peak_search_start, peak_search_end]
				function_fit, number of function to use to fit data in
				0: monoexponential, 3:biexponential, """
	
	set_params(params); 
	stf.measure();
	#this is in samples
	peak_index = stf.peak_index();
	stf.set_fit_start(peak_index, is_time=False);
	stf.set_fit_end(fit_end_time, is_time=True);
	fit_func = stf.leastsq(function_fit);
	return(fit_func)  
	
def jjm_calc_baseline_from_fit(fit_func, time_from_peak):

	#generate fit equation from parameters 
	time_base = np.arange(0, 40, stf.get_sampling_interval())
	
	decay = ((fit_func['Amp_0']+fit_func['Amp_1'])*(math.e**-(((1/fit_func['Tau_0'])+(1/fit_func['Tau_1']))*time_base)))+fit_func['Offset']
	
	baseline = decay[int(time_from_peak/stf.get_sampling_interval())]
	
	return(baseline)
	

				
def jjm_EPSC_peak_2(params_in, fit_end_time, function_fit, time_from_peak):
	"""params = [base_s, base_e, e_1_s, e_1_f, e_2_s, e_2_f]"""
	
	si = stf.get_sampling_interval()
		
	EPSC_1_peak = jjm_peak(params_in[0], params_in[1], params_in[2], params_in[3])
	
	fits = jjm_fit([params_in[0], params_in[1], params_in[2], params_in[3]], fit_end_time, function_fit)
	
	EPSC_2_baseline = jjm_calc_baseline_from_fit(fits, time_from_peak)
	
	peak_slice_start = int(params_in[4]/si)
	peak_slice_end = int(params_in[5]/si)
	
	EPSC_2_peak = np.min(stf.get_trace()[peak_slice_start:peak_slice_end])
	
	EPSC_rel = abs(EPSC_2_peak-EPSC_2_baseline)
	
	return(EPSC_1_peak, EPSC_rel)
	

def pp40_from_fits(params_in_i, fit_end_time_i, function_fit_i, time_from_peak_i):	
	
	EPSC_1 = []
	EPSC_2 = []
	pp40 = []
	for sweep in range(stf.get_size_channel()):
		stf.set_trace(sweep)
		output = jjm_EPSC_peak_2(params_in_i, fit_end_time_i, function_fit_i, time_from_peak_i)
		EPSC_1.append(output[0])
		EPSC_2.append(output[1])
		pp40.append(output[1]/output[0])
	df_out = pd.DataFrame({'EPSC_1':EPSC_1, 'EPSC_2':EPSC_2, 'pp40': pp40})
	
	fname = stf.get_filename().rstrip('abf')
	df_out.to_excel((fname+'pp40_from_fits.xlsx'))
		
	return(df_out)
		
def jjm_capacitance(params, function_fit):
	"""returns whole cell capacitance estimate (in pF)
		Inputs: params, list of [baseline_start, baseline_end, peak_search_start, peak_search_end]
				function_fit, number of function to use to fit data in
				0: monoexponential, 3:biexponential, """
	
	##call jjm_fit to get function parameters
	##jjm_fit returns a dictionary of fit parameters
	##{'SSE': 0.8533321749496918, 'Tau_0': 1.1497377662822952, 'Amp_0': 904.9731610575948, 'Offset': -275.5736999511719}
	##from AMP and tau calc capacitance
	
	return() 
	
def reimman_sum(boundary_params, baseline_value):
	"""calculates a reimann sum by multiplying the sampling interval times the amplitude of trace
	at each point, minus baseline
	Inputs: boundary_params(list of start and stop for sum), value for baseline
	"""
	reimann_start = boundary_params[0]/stf.get_sampling_interval();
	reimann_end = boundary_params[1]/stf.get_sampling_interval();
	integral_region = stf.get_trace()[int(reimann_start):int(reimann_end)];	
	integral_region_minusbaseline = integral_region-baseline_value;
	#y value is sampling interval
	delta_t = stf.get_sampling_interval(); 
	reimann = delta_t*integral_region_minusbaseline; 
	reimann_sum = sum(reimann); 	
	return(reimann_sum)
	
def capacitance_reimann(params, v_step):
	"""
	Inputs: params - [baseline_search_start, baseline_search_end, peak_search_start, peak_search_end]
	"""
	set_params(params); 
	stf.measure(); 
	baseline = stf.get_base(); 
	#peak index returns value in samples
	peak_index = stf.peak_index();
	peak_time = peak_index*stf.get_sampling_interval(); 
	reimann_params = [peak_time, params[3]];
	
	#reimann sum gives estimate for area under the curve, "charge accumulating on membrane capacitor"
	Q = reimman_sum(reimann_params, stf.get_base()); 
	#divide by delta V to get capacitance
	c = Q/v_step
	
		
	return(c)
	
def input_resistance(ss_end,v_amp):
	"""calculates input resistance using steady state current and voltage amplitude pulse Ri = (Voltage pulse in mV)/((SS_current - baseline) in pA)
	Inputs: ss_end (time to end search for steady state current), voltage amplitude
	Outputs: (SS_currrent - baseline), calculated Ri 
	"""
	
	si = stf.get_sampling_interval();
	#period to average for current is a 2ms region at end of pulse
	ss_start_samples = int((ss_end-2)/si) ; 
	ss_end_samples = int((ss_end)/si) ; 
	
	#gets steady state current values
	SS_absolute = np.mean(stf.get_trace()[ss_start_samples:ss_end_samples]);
	SS = SS_absolute-stf.get_base(); 
	
	#calculates input resistance
	Ri = v_amp/SS ;
	
	return(SS, Ri)
	
def add_capacitance_ri(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude):

	"""inputs: (baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude, EPSC1_s, EPSC1_e, EPSC2_s, EPSC2_e, sweep_start, sweep_end)
	output: numpy array where 1st column is capacitance transient amplitude, 2nd is series resistance, 3rd is 1st EPSC, 4th is 2nd EPCSC
	also writes output to .csv file"""
	
	## TODO: alter to have option for analyzing input number of sweeps
	num_sweeps = stf.get_size_channel(); 
	#print('there are')
	#print(num_sweeps)
	#print('sweeps in recording')
	#print('analyzing sweeps')
	#print(sweep_start)
	#print('to')
	#print(sweep_end)
	#sweeps_to_analyze = sweep_end - sweep_start
	
	#create array for results
	data_array = np.zeros((num_sweeps, 5)) ; 

	y = 0 
	for x in range(0, num_sweeps):
		#moves to next trace
		stf.set_trace(x);
		
		#jjm_resistance to get series resistance values 
		[cap_trans_amplitude, series_resistance] = jjm_resistance(baseline_start, baseline_end, cap_trans_start, cap_trans_end, amplitude);
		data_array[y][0] = cap_trans_amplitude;
		data_array[y][1] = series_resistance; 
		
		#input resistance calculation
		[SS, Ri] = input_resistance(cap_trans_end, amplitude);
		data_array[y][2] = SS; 
		data_array[y][3] = Ri; 
		
		#jjm_capacitance to get whole cell capacitance estimate
		whole_cell_capacitance = capacitance_reimann([baseline_start, baseline_end, cap_trans_start, cap_trans_end], amplitude);
		data_array[y][4] = whole_cell_capacitance; 
		y = y+1; 
			
	#print first few entries to check accuracy
	print(data_array[:3]); 
	
	#make data frame with array
	column_names = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)']; 
	data_array_df = pd.DataFrame(data_array, columns=column_names);
	
	#make csv file with data
	#file_name = stf.get_filename(); 
	#data_array_df.to_csv(file_name + '_stimfitanalysis.csv', sep=',')
		
	return(data_array_df)

def batch_add(directory, files, params):
	
	for fname in files:
		file_to_load = directory + fname[:12];
		print(file_to_load);
		stf.file_open(file_to_load);   
		#standard params = [100, 200, 0, 45, 10]

		cap_df = add_capacitance_ri(params[0], params[1], params[2], params[3], params[4]); 
		column_names = ['CapTransPeak(pA)', 'Rs(Mohm)', 'EPSC1', 'EPSC2']; 
		orig_df = pd.read_csv(file_to_load + '_stimfitanalysis.csv', names = column_names); 
		#pull out just EPSC values 
		to_combine = orig_df.ix[:,'EPSC1':'EPSC2'];

		#check if shapes not equal
		if cap_df.shape[0] != to_combine.shape[0]: 
			new_index = pd.Series(np.arange(cap_df.shape[0]-to_combine.shape[0], cap_df.shape[0]), name='new index'); 
			to_combine_ = pd.concat([to_combine, new_index], axis=1);
			to_combine = to_combine_.set_index('new index');
			
		combined_df = pd.concat([cap_df, to_combine], axis=1);
		
		##now resave combined as .csv
		combined_df.to_csv(directory+fname, sep=',');
	
	return(combined_df)
	
def open_file_list(file_list):
	for f_name in file_list:
		stf.file_open(f_name)
	return()
		
	
	
	
	
	
	
	
	
	
	
	
	
	