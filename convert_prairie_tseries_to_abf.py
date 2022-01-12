import csv
import numpy as np 
import os
import glob
import sys
import fnmatch 

"function for plotting .csv recording file from prairie view" 
"when calling script from command line, filename is 1st argument"

#code for finding values in .xml file
def get_value_from_line(line):
	value_list = []
	for char in line: 
		if (char.isdigit() == True) or (char == '.'):
			value_list.append(char) ; 
	value = ''.join(value_list); 
	return(value)

#code for loading channels of csv file and adjusting values based on gain/multipliers in xml file

def read_and_adjust_csv_recording(csv_filename, plot_var):
	print(csv_filename); 

# creates lists for time index and channel_0 values
	channel_0 = [] ;
	time = [] ; 

# imports channel_0 to list
	with open (csv_filename, 'rb') as csvfile:
		wholefile = csv.reader(csvfile, delimiter=',', quotechar='|')
	
		for row in wholefile:	
			time.append(row[0]);
			channel_0.append(row[1]) ;

	# need to find associated xml file 
	filename = csv_filename.rstrip('.csv'); 
	xml_file = filename + '.xml'
	print('xml_file is')
	print(xml_file)

	f = open(xml_file, 'r') ; 

	with open(xml_file) as f:
		xml_lines = [x.strip('\n') for x in f.readlines()] ; 

	for line in range(len(xml_lines)):
		## finds section for Input 0
		if '<Name>Input 0' in xml_lines[line]:	
		## divisor is located 9 lines up
		## gain is located 5 lines up 
			print 'reading info for channel:', xml_lines[line][14:21];
			gain = float(get_value_from_line(xml_lines[line+2])); 
			print 'gain is:', gain ;
			divisor = float(get_value_from_line(xml_lines[line+6]))
			print 'divisor is:', divisor;
			print(xml_lines[line]);
			print(xml_lines[line + 2]) ;
			print(xml_lines[line + 6])
			#channel_info =  
	

	# adjusts channel_0 per divisor in .xml file
	for i in range(len(channel_0)-1):
		channel_0[i+1] = float(channel_0[i+1]) ;
		channel_0[i+1] = (channel_0[i+1])/divisor ; 

	#prints 1st 10 entries
	print(channel_0[0:10]);

	# converts to floats
	for row in range(1,len(channel_0)):
		channel_0[row] = float(channel_0[row]);
		time[row] = float(time[row]);

	 

	sweep_array = np.array([time[1:], channel_0[1:]])
	#this should be a numpy array 
	# creates additional csv file for analysis
	newfilename =  filename + '_adjustedforimport.csv' ; 

	writer = csv.writer(open(newfilename, 'wb')) ;
	for i in range(len(channel_0)):
		writer.writerow([time[i], channel_0[i]]) ;
	return(sweep_array)
