import numpy as np
import stf

def remove_artifacts(first_artifact_start, length_time_to_remove, time_between_artifacts, number_of_artifacts):
	
	#get sampling interval of current tract
	
	sampling_interval = stf.get_sampling_interval()
	
	#convert input paramenters (units of time) to samples
	
	first_artifact_start_samples = float(first_artifact_start)/float(sampling_interval)
	
	length_time_to_remove_samples = float(length_time_to_remove)/float(sampling_interval)
	
	time_between_artifacts_samples = float(time_between_artifacts)/float(sampling_interval)
	
	#create variable for artifact end 
	
	first_artifact_end_samples = first_artifact_start_samples+length_time_to_remove_samples
	
	#get trace, convert trace to normal list
	
	original_trace = list(stf.get_trace())

	#create list for section of trace up until 1st artifact
	
	trace_artifacts_removed = original_trace[0:int(first_artifact_start_samples)]
	
	#add remaining sections of trace with artifacts removed
	
	for artifact_number in range(0, number_of_artifacts):
		
		start_sample = int(first_artifact_end_samples+(time_between_artifacts_samples*artifact_number))
		end_sample = int(first_artifact_start_samples+(time_between_artifacts_samples*(artifact_number+1)))
		
		print(start_sample)
		print(end_sample)
		
		trace_artifacts_removed.extend(original_trace[start_sample:end_sample])
		
	#plots trace in new window
	stf.new_window(trace_artifacts_removed)
	
	return(trace_artifacts_removed)
	
def remove_artifacts_from_sweeps(artifact_start_time, artifact_end_time):
	
	sampling_interval = stf.get_sampling_interval()
	artifact_start = int(artifact_start_time/sampling_interval)
	artifact_end = int(artifact_end_time/sampling_interval)
	
	continuous_trace = []
	output_artifacts_removed = []
	
	for sweep in range(stf.get_size_channel()):
		sweep_trace_before_artifact = stf.get_trace(sweep)[0:artifact_start]
		sweep_trace_after_artifact = stf.get_trace(sweep)[artifact_end:]
		sweep_trace = np.append(sweep_trace_before_artifact, sweep_trace_after_artifact)
		output_artifacts_removed.append(sweep_trace)
		continuous_trace.extend(sweep_trace)
		
	stf.new_window_list(output_artifacts_removed) 
	
	return(continuous_trace)
	
def save_for_import(file_name):
	np.savetxt(str(file_name)+'.csv', stf.get_trace())
	return()
	