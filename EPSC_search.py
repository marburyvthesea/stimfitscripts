import sys
add_to_path = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/wx-3.0-osx_cocoa', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks/stimfit', '.', '', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Applications/MacPorts/stimfit.app/Contents/Frameworks/stimfit', '/Applications/MacPorts/stimfit.app/Contents/Frameworks', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Library/Python/2.7/site-packages']
more_add_to_path = ['', '/Library/Python/2.7/site-packages/openpyxl-2.3.0-py2.7.egg', '/Library/Python/2.7/site-packages/et_xmlfile-1.0.1-py2.7.egg', '/Library/Python/2.7/site-packages/jdcal-1.2-py2.7.egg', '/Library/Python/2.7/site-packages/XlsxWriter-0.7.7-py2.7.egg', '/Library/Python/2.7/site-packages', '/Library/Python/2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC', '/Library/Python/2.7/site-packages'];
for path in add_to_path:
	sys.path.append(path);  
for path in more_add_to_path:
	sys.path.append(path);
import numpy as np
import stf
import math
import pandas as pd
import scipy.stats

class event_counts(object):

##create class for event count file 
#make this a container to store the name of the saved event files, the whole trace, 
#functions will be options to calculate properties of detected EPSCs 
#with option to output spreadsheet file
	def __init__(self, event_file_name, whole_trace_file_name, sweep):
		self.event_file = event_file_name
		self.whole_trace_file = whole_trace_file_name
		self.sweep_number = sweep
		
	def get_time_points(self):
	#need to get specified sweep of the whole trace to an array
		stf.file_open(self.whole_trace_file)
		stf.set_trace(self.sweep_number)
		trace_array = stf.get_trace()
	
	#runs program to find sample indicies of selected EPSCs and converts to times
		stf.file_open(self.event_file)
		event_samples_list = find_sample_points_of_detected_events(trace_array)
		sampling_interval = stf.get_sampling_interval()
		event_times_list = [(sample*sampling_interval) for sample in event_samples_list]
		
		return(event_times_list)
		
	def output_EPSC_info(self, filename_for_output):
		#runs get times
		times = self.get_time_points()
		#creates data frame for writing to excel
		mEPSC_data = pd.DataFrame({'EPSC times':times})
		mEPSC_data.to_excel(filename_for_output)
		#np.savetxt(times, self.event_file+'event_times')
		return()
		
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
	
	while found == False:
		
		segment_to_test = full_trace[start_index:(start_index+trace_to_find_length)]
		
		
		if list(segment_to_test) == list(trace_to_find):
			found = True
			found_index = start_index
		
		
		start_index = start_index+1
				
		
	return(found_index)

def run_automated_search_triexponential(dictionary_inputs, correlation_coeff_threshold):
	
	trace_region_to_search = (0, stf.get_size_trace)
	search_period = 200
	threshold = correlation_coeff_threshold
	min_btw_events = 1
	tau_rise = dictionary_inputs['Tau_0']
	tau_1_decay = dictionary_inputs['Tau_1']
	tau_2_decay = dictionary_inputs['Tau_2']
	


def automated_search_triexponential(trace_region_to_search, search_period, threshold, min_btw_events, tau_rise, tau_1_decay, tau_2_decay):
	"""searches section of trace based on a user input triexponential function (tau_rise, tau_1_decay, tau_2_decay)"""
	#converts some inputs to sample points
	min_samples_btw_events = min_btw_events/stf.get_sampling_interval()
	
	#pull out region to search
	region_to_search = stf.get_trace()[trace_region_to_search[0]:trace_region_to_search[1]]
	
	#list to store detected events
	event_times = []
	
	#creates vector of time points 
	t = np.linspace(0,50,(50/stf.get_sampling_interval()))
	
	#creates triexponential pattern function 
	p_t = [(1-math.exp(-(t_point-0)/tau_rise))*(math.exp(-(t_point-0)/tau_1_decay))*(math.exp(-(t_point-0)/tau_2_decay)) for t_point in t]
	
	#slides window along 
	pt = 0 
	while pt < range(len(region_to_search)-int(min_samples_btw_events)):
		
		EPSC_test = stf.get_trace()[pt:(pt+(search_period/stf.get_sampling_interval()))]
		
		corr_coeff = stats.pearsonr(p_t, EPSC_test)[0]
		
		if corr_coeff > threshold:
			
			stf.set_marker(pt, region_to_search[trace_region_to_search[0]+pt])
			
			event_times.append(pt*stf.get_sampling_interval())
			
			pt += min_samples_btw_events
			
		else:
			
			pt += 1
			
	return(event_times)
			
def baseline_from_linear_regression():

	y_values_trace = stf.get_trace()
	x_values_trace = range(0, len(stf.get_trace()))
	results = scipy.stats.linregress(x_values_trace, y_values_trace)
	
	return(results)
		
def plot_linear_regression(slope, intercept, length):
	
	x_values = range(0, length, 1)
	y_values = [((slope*x) + round(intercept)) for x in x_values]
	
	regression_plot = np.array([x_values, y_values])
	
	stf.new_window_matrix(regression_plot)
	
	return(x_values, y_values)
	
def run_plot_linear_regression():
	lin_results = baseline_from_linear_regression()
	plot_linear_regression(lin_results[0], lin_results[1], len(stf.get_trace()))
	
	return()
	
			
		