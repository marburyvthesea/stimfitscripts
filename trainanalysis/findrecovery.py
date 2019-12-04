


def find_trains(directory):
	## find recovery files matching peak files
	peak_files = [file for file in os.listdir(directory) 
	if '_peaks.csv' in file]
	
	#finds this stimfit analysis recovery files (whether or not they exist) in directory based on analyzed peak_files
	#in directory
	recovery_files = [file[:7]+str(int(file[7])+1)+'.abf_stimfitanalysis.csv' for file in peak_files]
	
	baseline_files = [file[:7]+str(int(file[7])-1)+'.abf_stimfitanalysis.csv' for file in peak_files]
	
	return(peak_files, recovery_files, baseline_files)