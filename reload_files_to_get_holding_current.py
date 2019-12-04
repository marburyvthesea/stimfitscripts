import sys
import os
add_to_path = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/wx-3.0-osx_cocoa', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks/stimfit', '.', '', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Applications/MacPorts/stimfit.app/Contents/Frameworks/stimfit', '/Applications/MacPorts/stimfit.app/Contents/Frameworks', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Library/Python/2.7/site-packages']
more_add_to_path = ['', '/Library/Python/2.7/site-packages/openpyxl-2.3.0-py2.7.egg', '/Library/Python/2.7/site-packages/et_xmlfile-1.0.1-py2.7.egg', '/Library/Python/2.7/site-packages/jdcal-1.2-py2.7.egg', '/Library/Python/2.7/site-packages/XlsxWriter-0.7.7-py2.7.egg', '/Library/Python/2.7/site-packages', '/Library/Python/2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC', '/Library/Python/2.7/site-packages'];
for path in add_to_path:
	sys.path.append(path);  
for path in more_add_to_path:
	sys.path.append(path);
import pandas as pd	
import pandas as pd
import numpy as np
import xlrd
import xlsxwriter
import sys
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/jupyterscripts/')
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts/')
import compile_by_abf_file as cf
import stf
import prairie_tseries_toabf_stimfit as pv

def search_dir_iterative_for_extension_tupleoutput(rootDir, extension_to_search):
    ##Idea is to compile a list of all my abf files that can be referenced against the stimfit analysis files
    """Function to search a directory and subdirectories for file names containing a given file extension 
    inputs: root directory, file extension(string at end of filename)
    output: list of tuples, filenames and associated directories"""
    # create a numpy array to append files and directories to
    tuple_array = np.array([('file name','directory')]); 
    for dirName, subdirList, fileList in os.walk(rootDir):
        #print('Found directory: %s' % dirName)
        files_with_string = [] ; 
        for fname in fileList:
            if fname.endswith(extension_to_search):
                files_with_string.append((fname, dirName)); 
        #print('found ' + str(len(files_with_string)) + ' files')            
        #only add to list if a file is found
        if len(files_with_string)>0:
            np_array = np.array(files_with_string);
            tuple_array = np.vstack((tuple_array, np_array)); 
    return(tuple_array);



def return_holding_current(selected_wb, raw_wb):
	
	wb_ = xlrd.open_workbook(selected_wb)
	sheets = [str(sheet.name) for sheet in wb_.sheets() if 'normalized' in str(sheet.name)]
	sheets_to_load_from_raw = [sheet.strip('normalized')+'iled_files' for sheet in sheets]
	wb_raw = xlrd.open_workbook(raw_wb)
	raw_sheets = [str(sheet.name) for sheet in wb_raw.sheets() if 'compiled_files' in str(sheet.name)]
	baseline_files = []
	exp_files = []
	for raw_sheet in raw_sheets:
		
		try:
			df = pd.read_excel(wb_raw, engine='xlrd', sheetname=str(raw_sheet), index_col = [0,1])
			
			
			files_in_experiment = df.index.levels[0].values[:2]
			
			baseline_files.append(files_in_experiment[0])
			exp_files.append(files_in_experiment[1])
		except:
			pass
	
	print(baseline_files)
	
	baseline_files_to_load = []
	exp_files_to_load = []
	rootsearchdir = '/Users/johnmarshall/Documents/Analysis/RecordingData/'
	for fname, efname in zip(baseline_files, exp_files) : 
		try:
			if 'TSeries' in fname:
				dir_data = search_dir_iterative_for_extension_tupleoutput(rootsearchdir, str(fname+'vrecd_loaded.csv'))
			else:
				dir_data = search_dir_iterative_for_extension_tupleoutput(rootsearchdir, fname)
				baseline_files_to_load.append(dir_data[1][1]+'/'+dir_data[1][0])
				dir_data = search_dir_iterative_for_extension_tupleoutput(rootsearchdir, efname)
				exp_files_to_load.append(dir_data[1][1]+'/'+dir_data[1][0])
		except IndexError:
			print('could not find:', fname)
		pass
	
	print(baseline_files_to_load)
	dict_to_return = {}
	holding_current_time_series = {}
	for f, ef in zip(baseline_files_to_load, exp_files_to_load):
		currents = []
		print(f)
		
		if f.endswith('.abf'):
			stf.file_open(f)
		else:
			if 'TSeries' in f:
				sweeps_compiled_from_pv_tseries = pv.import_t_series_episodic(f.strip('vrecd_loaded.csv'))
				pv.plot_episodic_array(sweeps_compiled_from_pv_tseries) 
		
		baselines = []
		for sweep in range((stf.get_size_channel()-30),stf.get_size_channel()):
			stf.set_trace(sweep)
			stf.set_base_start(100, is_time=True)
			stf.set_base_end(100, is_time=True)
			baselines.append(stf.get_base())
		file_baseline = np.mean(baselines)
		currents.append(file_baseline)
		holding_current_time_series[f] = baselines
		
		stf.file_open(ef)
		baselines = []
		for sweep in range((stf.get_size_channel()-15),stf.get_size_channel()):
			stf.set_trace(sweep)
			stf.set_base_start(100, is_time=True)
			stf.set_base_end(100, is_time=True)
			baselines.append(stf.get_base())
		file_baseline = np.mean(baselines)
		currents.append(file_baseline)
		dict_to_return[f] = currents
	
	df_out = pd.DataFrame(dict_to_return)
	time_series_df = pd.DataFrame(holding_current_time_series)
	
	df_out.to_excel('/Users/johnmarshall/Documents/Analysis/eCB_paper/holding_current.xlsx')
	time_series_df.to_excel('/Users/johnmarshall/Documents/Analysis/eCB_paper/holding_current_time_series.xlsx')
	
	return(df_out)

def return_base_for_file(start_sweep, end_sweep):
	
	#dict_to_return = {}
	baselines = []
	for sweep in range(start_sweep, end_sweep):
		stf.set_trace(sweep)
		stf.set_base_start(100, is_time=True)
		stf.set_base_end(125, is_time=True)
		baselines.append(stf.get_base())
	file_baseline = np.mean(baselines)
	#dict_to_return[stf.get_filename()] = file_baseline
	
	#df_out = pd.DataFrame(dict_to_return)
	#file_name = stf.get_filename()
	#df_out.to_excel('/Users/johnmarshall/Documents/Analysis/eCB_paper/'+str(file_name)+'holding_current.xlsx')
	
	return(file_baseline)
	






