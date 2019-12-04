import sys
import numpy as np
import stf
import spells


def find_AP_peaks(start_msec, delta_msec, current_start, current_delta, threshold_value, deflection_direction, mark_option):
	
	""" count number of APs in traces with current injection/gradually increasing steps
	inputs: (time (msec) to start search, length of search region, starting current value, current delta between traces, threshold value, deflection direction ('up'/'down'), mark traces (True/False))"""	
	
	event_counts = np.zeros((stf.get_size_channel(),2)); 
	
	for trace_ in range(stf.get_size_channel()):
		event_counts[trace_][1] = spells.count_events(start_msec, delta_msec, threshold=threshold_value, up=deflection_direction, trace=trace_, mark=mark_option); 
		event_counts[trace_][0] = current_start + (current_delta*trace_) ; 
		
	loaded_file = stf.get_filename()[:-4] ; 
	np.savetxt(loaded_file + '_AP_counts.csv', event_counts, delimiter=',', newline='\n'); 
	return(event_counts)
	
def find_AP_peak_ADP(start_msec, delta_msec, current_start, current_delta, threshold_value, deflection_direction, mark_option):
	
	""" count number of APs in traces with current injection/gradually increasing steps
	inputs: (time (msec) to start search, length of search region, starting current value, current delta between traces, threshold value, deflection direction ('up'/'down'), mark traces (True/False))"""	
	
	loaded_file = stf.get_filename()[:-4] 
	event_counts = np.zeros((stf.get_size_channel(),2)); 	
	for trace_ in range(stf.get_size_channel()):
		##gets AP counts and sample points in current trace 
		if deflection_direction == 'up':
			direction_input = True
		else:
			direction_input = False
		[trace_count, sample_points, time_points] = jjm_count(start_msec, delta_msec, threshold=threshold_value, up=direction_input, trace=trace_, mark=mark_option); 
		
		##gets ADP values 
		values, indicies = find_ADPs(sample_points)
		print(values)
		out_array = np.array([values, indicies])
		
		np.savetxt(loaded_file + 'trace' + str(str(trace_).zfill(3)) +'ADP_values.csv', out_array, delimiter=',', newline='\n')
		
		event_counts[trace_][1] = trace_count
		
		event_counts[trace_][0] = current_start + (current_delta*trace_) ; 
		
	np.savetxt(loaded_file + '_AP_counts.csv', event_counts, delimiter=',', newline='\n'); 
	return(event_counts)

	
	
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
	count_events(500,1000,0,False,-10,i) returns the number of events found below -10 in the 
	trace i and shows the corresponding stf markers. """
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
	
	time_points = [sample_point*dt for sample_point in sample_points];
	return (EventCounter, sample_points, time_points)
	
def find_ADPs(AP_peak_indicies):
	ADP_values = []
	ADP_indicies = []
	##slices 
	for peak in range(len(AP_peak_indicies)-1):
		ADP_search = stf.get_trace()[AP_peak_indicies[peak]:AP_peak_indicies[peak+1]]
		min_value = np.min(ADP_search)
		min_index = AP_peak_indicies[peak] + np.argmin(ADP_search)
		stf.set_marker(min_index, min_value)
		ADP_values.append(min_value)
		ADP_indicies.append(min_index)
			
	return(ADP_values, ADP_indicies)













