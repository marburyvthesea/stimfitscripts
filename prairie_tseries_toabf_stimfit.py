import os
import glob
import sys 
import convert_prairie_tseries_to_abf
import numpy as np 
import stf

"function will find the newest t-series file, compile a list of the voltage_recording .csv files"
"and then call read csv to convert, followed by find peaks, 1st tuple passed is times for baseline, followed by variable for peaks"
"it should output the peaks to one central csv file" 
"has functions to load and plot into stimfit"
#ask if plot
plot = 'n'
if plot == 'y':
	plot_var = True
else:
	plot_var = False

def get_sampling_delta(time_vector):
	print(time_vector[1])
	print(time_vector[0])
	sampling_delta = float(time_vector[1]-time_vector[0]) ; 
	return(sampling_delta)
	
def get_samples_from_time(times_input, sampling_delta):
	times = [] 
	for x, y in zip(range(0, len(times_input), 2),range(1, len(times_input), 2)):
		time_tuple = tuple((float(times_input[x])/sampling_delta, float(times_input[y])/sampling_delta)); 
		times.append(time_tuple) ; 
	return(times)	

## function to find csv files in directory
def compile_csv_files(directory):
	os.chdir(directory)
	csv_file_list = [];
	for file in os.listdir(os.curdir):
		if file.endswith('VoltageRecording_001.csv'):
			csv_file_list.append(file);
	
	return(csv_file_list);


## Create arrays for whole sweep
def import_t_series_episodic(input_directory):
	
	file_list = compile_csv_files(input_directory);
	
	#gets time from 1st file 
	episodic_sweeps_array = convert_prairie_tseries_to_abf.read_and_adjust_csv_recording(file_list[0], plot_var)[0] ; 

	#makes array with entire sweeps appended to time vector 
	for file in file_list:
		episodic_sweep = convert_prairie_tseries_to_abf.read_and_adjust_csv_recording(file, plot_var)[1]; 
		
		episodic_sweeps_array = np.vstack((episodic_sweeps_array, episodic_sweep)); 
	
	#returns to previous directory
	os.chdir('..')
	
	np.savetxt(input_directory + 'vrecd_loaded.csv', np.transpose(episodic_sweeps_array), delimiter=',', newline='\n');

	plot_episodic_array(episodic_sweeps_array); 
	
	return(episodic_sweeps_array)
	
def plot_episodic_array(sweeps_array):
	#extracts time
	time = sweeps_array[0] ; 
	
	#gets sampling interval 
	sampling_interval = get_sampling_delta(time) ; 
	
	#plots all traces in single window
	stf.new_window_matrix(sweeps_array[1:]) ; 
	
	#adjusts sampling interval  
	stf.set_sampling_interval(sampling_interval); 
		
	return()