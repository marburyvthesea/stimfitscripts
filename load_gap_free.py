import sys
sys.path.append('/Users/johnmarshall/anaconda3/envs/stimfit_env/lib/python2.7/site-packages/')
from neo import io
import csv
import numpy as np 


def load_gap_free_trace(file_to_load):
	"""imports abf file, must be in directory of file
	Input: abf_file_to_load"""
	
	filename = file_to_load; 
	experiment_name = filename.rstrip('.abf');

	r = io.AxonIO(filename=file_to_load)
	#bl = r.read_block(lazy=False, cascade=True)
	bl = r.read_block(lazy=False)
	#segments are sweeps

	
	print bl.segments[0].analogsignals[0].magnitude
	
	##get sampling rate
	sampling_rate = bl.segments[0].analogsignals[0].sampling_rate
	print(sampling_rate)

	##adds channel 0 from each sweep to array 
	print('file has')
	print(len(bl.segments))
	print('sweeps')
	print(len(bl.segments[0].analogsignals[0].magnitude))
	print('samples')
	channel_array = np.empty((len(bl.segments)+1,(len(bl.segments[0].analogsignals[0])))); 
	print(channel_array.shape)
	for sweep in range(len(bl.segments)):
		channel_0_sweep = [] 
		for data_point in range(len(bl.segments[sweep].analogsignals[0].magnitude)):	
			#print(bl.segments[sweep].analogsignals[0].magnitude[data_point])
			channel_array[sweep+1][data_point] = (bl.segments[sweep].analogsignals[0].magnitude[data_point]);
	
	
	print channel_array[0][0:10]
	


	## make additional row for time
	samplingrate_Hz = sampling_rate.magnitude ;
	sampling_interval_msec = (1000 / float(samplingrate_Hz));
	for time_point in range(len(bl.segments[sweep].analogsignals[0].magnitude)):
		channel_array[0][time_point] = (float(time_point)*sampling_interval_msec);  

	## write a csv file 

	np.savetxt(experiment_name + 'abf_to_csv.csv', np.transpose(channel_array), delimiter=',', newline='\n');
	return(channel_array)


 



