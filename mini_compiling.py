import sys
sys.path.append('/Library/Python/2.7/site-packages');
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
sys.path.append('/anaconda3/envs/Py27/lib/python2.7/site-packages')
from scipy import ndimage
import stf
import numpy as np
import pandas as pd
import os

##cumulative distribtuions
#import scipy.stats as stats
#import matplotlib.pyplot as plt
#import numpy as np

#x = np.linspace(0,20,100)
#cdf = stats.binom.cdf
#plt.plot(x,cdf(x, 50, 0.2))
#plt.show()



def get_extracted_event_files(bottom_directory):
    files_to_return = []
    for root, dirs, files in os.walk(bottom_directory, topdown=False):
        
        for name in files:
            if ('extracted' in name) and (name.endswith('.h5')):
                files_to_return.append(str(os.path.join(root, name)))
    return(files_to_return)

def batch_compile(file_list, summary_file_name):

	means_dict = {}
	
	for fname in file_list:
		stf.file_open(fname)
		
		file_df = compile_amplitudes_in_trace()
		
		means_dict[stf.get_filename()] = file_df.mean(axis=0)
		
	summary_df = pd.DataFrame(means_dict)
	summary_df.to_excel(str(summary_file_name) + '_amplitudes_compiled.xlsx')
	
	return(summary_df)
		

def compile_amplitudes_in_trace():
	
	# for each trace in file run find_baseline_amplitudes
	output_array = np.array(['baseline', 'peak', 'peak_from_baseline'])
	
	for trace in range(stf.get_size_channel()):
		
		stf.set_trace(trace)
		
		fba_output = find_baseline_amplitude(10)
		
		output_array = np.vstack([output_array, fba_output])
		
	output_df = pd.DataFrame(output_array[1:], columns = output_array[0], dtype=float)
		
	output_df.to_excel(str(stf.get_filename().rstrip('.h5')+'_amplitudes.xlsx'))
	
	return(output_df)		

def find_baseline_amplitude(sigma):
	# gaussian filter with sigma 10
	trace_ = stf.get_trace()
	trace_filtered = ndimage.filters.gaussian_filter(trace_, sigma)
	# take derivative
	si = stf.get_sampling_interval()
	#read V values from trace, 
	V_values = stf.get_trace()
	#compute dv and by iterating over voltage vectors
	dv = [V_values[i+1]-V_values[i] for i in range(len(V_values)-1)]
	#compute dv/dt
	dv_dt = [(dv[i]/si) for i in range(len(dv))]
	# find index of derivative peak
	deriv_max = np.argmin(dv_dt)
	# use derivative peak index to get baseline from original trace
	# use a mean of 10 sample points
	baseline = np.mean(trace_[deriv_max-10:deriv_max])
	stf.set_marker(deriv_max, baseline)
	peak_amplitude = np.min(stf.get_trace())
	peak_index = np.argmin(stf.get_trace())
	stf.set_marker(peak_index, peak_amplitude)
	peak_from_baseline = peak_amplitude-baseline
		
	return(baseline, peak_amplitude, peak_from_baseline)

def write_dictionaries(full_file):
	dicts_list = []
	#find current directory
	current_dir = os.getcwd() 
	for fname in os.listdir(current_dir):
		if (full_file[:8] in fname) and (fname.endswith('h5')):
			dicts_list.append({full_file : [fname]})
	
	return(dicts_list)
		
def batch_find_ISIs(dictionary_for_traces):
# instead of dictionaries for traces have everything loaded into an excel or csv file, first column contains the
# name of the full trace, second column contains the name of the extracted event files
# then load each line of the column and run this function on the dictionary 
# so modify this to save the df based on extracted file name 
	
	output_dfs = []
	whole_trace_file_names = []
	
	for whole_trace_file_name in dictionary_for_traces.keys():
		
		event_file_dictionary_ISI = {} 
		
		for event_file_name in dictionary_for_traces[whole_trace_file_name]:
			
			sweep_num = event_file_name[16]
			
			print('sweep' + sweep_num)
					
			ISI_index = find_sample_points_of_detected_events(whole_trace_file_name, event_file_name, int(sweep_num))
			
			event_file_dictionary_ISI[event_file_name] = ISI_index
		
		event_file_df = pd.DataFrame(event_file_dictionary_ISI)
				
		output_dfs.append(event_file_df)
		whole_trace_file_names.append(whole_trace_file_name)
	
	fname_to_save = dictionary_for_traces.values[0]
		
	output_df = pd.concat(output_dfs, keys = whole_trace_file_names)
	output_df.to_excel(fname_to_save)
								
	return(output_df)
			
def find_sample_points_of_detected_events(whole_trace_file, extracted_events_file, sweep_num):
	"""takes the window of detected events from stimfit and, for each events, runs through the full trace to pull out time (in samples) of event
	"""
	#open and load trace from whole file
	stf.file_open(whole_trace_file)
	stf.set_trace(sweep_num)
	whole_trace = stf.get_trace()
	sampling_interval = stf.get_sampling_interval()
	
	#open extracted events file
	stf.file_open(extracted_events_file)
	
	time_points = []
	for trace in range(stf.get_size_channel()):
		stf.set_trace(trace)
		trace_to_search = stf.get_trace(trace)
		# run find trace with updated search index
		# start at sample = 0 for first run through
		if len(time_points) == 0:
			sample_start = 0
		else:
			sample_start = int(time_points[len(time_points)-1]/sampling_interval)
			
		output_index = sub_func_find_trace(trace_to_search, whole_trace, sample_start)
		time_point = output_index*sampling_interval
		time_points.append(time_point)
		
	return(time_points)

def sub_func_find_trace(trace_to_find, full_trace, start_sample):
	"""finds start point of single trace in full trace
	inputs: trace_to_find, full_trace as numpy arrays, start_sample	for search
	NB: this runs slow, need to find improvements to searching whole trace
	
	"""	
	trace_to_find_length = trace_to_find.size
	full_trace_length = full_trace.size
	
	start_indicies = []
	#print(trace_to_find)
	
	found = False
	start_index = start_sample
	
	found_index = 'event not found'
	
	while (found == False) and (start_index<=full_trace_length):
		
		segment_to_test = full_trace[int(start_index):int((start_index+trace_to_find_length))]
		
		
		if list(segment_to_test) == list(trace_to_find):
			found = True
			found_index = start_index
		
		
		start_index = start_index+1
				
	
		
	return(found_index)

	
	
	

