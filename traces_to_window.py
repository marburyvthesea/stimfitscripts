import stf

def get_traces(start, end):
	trace_list = []; 
	for x in range(start, end):
		trace = stf.get_trace(x); 
		new_window(trace); 
		
	return()