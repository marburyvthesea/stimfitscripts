import os 


def find_trains_iterative(root_directory):

	train_files = []
	
	for root, dirs, files in os.walk(root_directory, topdown=False):
		
		for fname in files:
	
			if '_peaks.csv' in fname:
				train_files.append([root , fname])
	
	return(train_files)

def find_trains(directory):
	## find recovery files matching peak files
	peak_files = [file for file in os.listdir(directory) 
	if '_peaks' in file]
	
	#finds this stimfit analysis recovery files (whether or not they exist) in directory based on analyzed peak_files
	#in directory
	recovery_files = [file[:7]+str(int(file[7])+1)+'.abf_stimfitanalysis.csv' for file in peak_files]
	
	baseline_files = [file[:7]+str(int(file[7])-1)+'.abf_stimfitanalysis.csv' for file in peak_files]
	
	return(peak_files, recovery_files, baseline_files)
	
def find_baseline_and_recovery_for_file(train_file, directory, searchdir):
	if train_file[5:8]>=1:
		recovery = train_file[:5]+str(int(train_file[5:8])+1).zfill(3)+'.abf_stimfitanalysis.csv'
		baseline = train_file[:5]+str(int(train_file[5:8])-1).zfill(3)+'.abf_stimfitanalysis.csv'
	for dirName, subdirList, fileList in os.walk(searchdir):
		if recovery in fileList:
			recovery_out = recovery
			recovery_dir = str(dirName)
	for dirName, subdirList, fileList in os.walk(searchdir):
		if baseline in fileList:
			baseline_out = baseline 
			baseline_dir = (dirName)
	
	if 'recovery_out' not in locals():
		recovery_out = 'did not find in' +str(searchdir)
	if 'recovery_dir' not in locals():
		recovery_dir = 'did not find in' +str(searchdir)
	if 'baseline_out' not in locals():
		baseline_out = 'did not find in' +str(searchdir)
	if 'baseline_dir' not in locals():
		baseline_dir = 'did not find in' +str(searchdir)
	
	return(recovery_out, recovery_dir, baseline_out, baseline_dir)

def create_compiled_df(train_file_df, file_type, column_name, file_list):
    """
    Inputs: 
    		train_file_df (df with structure: numeric indicies, 'trainfile' column with *peaks_normalized*.csv file,
    		'directory' column with location, 'recovery' coulmn with recovery file, 'baseline' column with baseline file,
    		'recovery_dir' and 'baseline' dir columns)
    		file_type ('recovery' or 'baseline')
            column_name to pull out (e.g. 'EPSC1' or 'Cap_TransPeak(pA)')
    Output: df with all applicable columns lined up"""
    ##TO DO: take a list of experiments as input to check against
    ##loop to pull out EPSCs and create data frame with all recovery EPSCs lined up
    time_series = pd.Series(np.arange(0,30,float(1)/float(3)))
    time_series_df = pd.DataFrame(time_series, columns=['time'])
    for x in train_file_df.index.values:
        if train_file_df.loc[x]['trainfile'][:8] in file_list:
        
            r_file = train_file_df.loc[x][file_type]
            if r_file.endswith('.csv'):
                to_load = train_file_df.loc[x][file_type+'_dir']+'/'+train_file_df.loc[x][file_type]
                recovery_df = pd.read_csv(to_load)
                #if file doesn't have column values as a header need to assign them
                if recovery_df.columns.values[0] != 'CapTransPeak(pA)':
                    recovery_df = pd.read_csv(to_load, names=['CapTransPeak(pA)', 'Rs(Mohm)', 'EPSC1', 'EPSC2'])
                #if file does not have columsn for SS, Ri, Capacitance need to inser those
                if recovery_df.columns.values[2] != 'SS(pA)':
                    recovery_df_full = recovery_df.reindex(index=recovery_df.index, columns = ['CapTransPeak(pA)','', 'Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)', 'EPSC1','', 'EPSC2'])
                else:
                    recovery_df_full = recovery_df
                EPSCs_to_concat = pd.DataFrame(recovery_df_full[column_name].as_matrix(), columns = [str(train_file_df.loc[x][file_type])])
                time_series_df = pd.concat([time_series_df, EPSCs_to_concat], axis=1)
    return(time_series_df)