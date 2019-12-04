import sys
sys.path.append('/Library/Python/2.7/site-packages');
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
from scipy import ndimage
import stf
import numpy as np
import pandas as pd
import os


def get_extracted_event_files():
	for root, dirs, files in os.walk(bottom_directory, topdown=False):
		files_to_return = []
		for name in files:
			if 'extraced' in name:
				print(os.path.join(root, name)) 
				files_to_return.append(name) 
	
	return(files_to_return)

def batch_compile(file_list, summary_file_name):

	means_dict = {}
	
	for fname in file_list:
		stf.file_open(fname)
		
		file_df = compile_amplitudes_in_trace()
		
		means_dict[stf.get_filename()] = file_df.mean(axis=0)
		
	summary_df = pd.DataFrame(means_dict)
	summary_df.to_excel(str(summary_file_name) + '.xlsx')
	
	return(summary_df)
		

def compile_amplitudes_in_trace():
	
	# for each trace in file run find_baseline_amplitudes
	output_array = np.array(['baseline', 'peak', 'peak_from_baseline'])
	
	for trace in range(stf.get_size_channel()):
		
		stf.set_trace(trace)
		
		fba_output = find_baseline_amplitude(10)
		
		output_array = np.vstack([output_array, fba_output])
		
	output_df = pd.DataFrame(output_array[1:], columns = output_array[0], dtype=float)
		
	output_df.to_excel(str(stf.get_filename()[-40:-3])+'.xlsx')
	
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
	peak_amplitude = np.min(stf.get_trace())
	peak_from_baseline = peak_amplitude-baseline
		
	return(baseline, peak_amplitude, peak_from_baseline)
	
def batch_find_ISIs(dictionary_for_traces):
	
	output_dictionary_ISI = {}
	
	for whole_trace_file in dictionary_for_traces.keys():
		
		event_file_dictionary_ISI = {} 
		
		for event_file in dictionary_for_traces[whole_trace_file]:
			
			stf.file_open(event_file)
		
			sweep_num = event_file[20]
			
			print('sweep' + sweep_num)
					
			ISI_index = find_sample_points_of_detected_events(whole_trace_file, int(sweep_num))
			
			event_file_dictionary_ISI[event_file] = ISI_index
			
		output_dictionary_ISI[whole_trace_file] = event_file_dictionary_ISI
							
	return(output_dictionary_ISI)
			
def find_sample_points_of_detected_events(whole_trace):
	"""takes the window of detected events from stimfit and, for each events, runs through the full trace to pull out time (in samples) of event
	"""
	trace_indicies = []
	for trace in range(stf.get_size_channel()):
		trace_to_search = stf.get_trace(trace)
		# run find trace with updated search index
		# start at sample = 0 for first run through
		if len(trace_indicies) == 0:
			sample_start = 0
		else:
			sample_start = trace_indicies[len(trace_indicies)-1]
			
		output_index = sub_func_find_trace(trace_to_search, whole_trace, sample_start)
		trace_indicies.append(output_index)
		
	return(trace_indicies)

def sub_func_find_trace(trace_to_find, full_trace_file, start_sample):
	"""finds start point of single trace in full trace
	inputs: trace_to_find, full_trace as numpy arrays, start_sample	for search
	NB: this runs slow, need to find improvements to searching whole trace
	
	"""	
	full_trace = full_trace_file

	trace_to_find_length = trace_to_find.size
	full_trace_length = full_trace.size
	
	start_indicies = []
	#print(trace_to_find)
	
	found = False
	start_index = start_sample
	
	while found == False:
		
		segment_to_test = full_trace[start_index:(start_index+trace_to_find_length)]
		
		
		if list(segment_to_test) == list(trace_to_find):
			found = True
			found_index = start_index
		
		
		start_index = start_index+1
				
		
	return(found_index)

	
	
	

