import sys
add_to_path = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/wx-3.0-osx_cocoa', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks/stimfit', '.', '', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Applications/MacPorts/stimfit.app/Contents/Frameworks/stimfit', '/Applications/MacPorts/stimfit.app/Contents/Frameworks', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Library/Python/2.7/site-packages']
more_add_to_path = ['', '/Library/Python/2.7/site-packages/openpyxl-2.3.0-py2.7.egg', '/Library/Python/2.7/site-packages/et_xmlfile-1.0.1-py2.7.egg', '/Library/Python/2.7/site-packages/jdcal-1.2-py2.7.egg', '/Library/Python/2.7/site-packages/XlsxWriter-0.7.7-py2.7.egg', '/Library/Python/2.7/site-packages', '/Library/Python/2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC', '/Library/Python/2.7/site-packages'];
for path in add_to_path:
	sys.path.append(path);  
for path in more_add_to_path:
	sys.path.append(path);
import jjm_analysis
import stf
import numpy as np
import pandas as pd

def fit_experiment(params, pulse_length, function_to_fit):
	
	num_sweeps = stf.get_size_channel(); 
	stf.set_channel(0); 
	stf.set_trace(0);
	
	#jjm_analysis.set_params(params); 
	#stf.measure();
	#this is in samples
	#peak_index = stf.peak_index();
	#stf.set_fit_start(peak_index, is_time=False);
	#fit_start_time = peak_index*stf.get_sampling_interval(); 
	#stf.set_fit_end(fit_start_time+pulse_length-(10*stf.get_sampling_interval()), is_time=True);
	#fit_func = stf.leastsq(function_to_fit);
	#fit_func['Baseline(pA)']=stf.get_base();
	#fit_df = pd.DataFrame(fit_func, index=[0]); 
		
	fits = [];
	traces = []; 
	for x in range(0, num_sweeps):
		stf.set_trace(x);
		jjm_analysis.set_params(params); 
		stf.measure();
		#this is in samples
		peak_index = stf.peak_index();
		stf.set_fit_start(peak_index, is_time=False);
		fit_start_time = peak_index*stf.get_sampling_interval(); 
		stf.set_fit_end(fit_start_time+pulse_length-(10*stf.get_sampling_interval()), is_time=True);
		sweep_fit = stf.leastsq(function_to_fit);
		sweep_fit['Baseline(pA)']=stf.get_base();
		fits.append(sweep_fit);
		traces.append(x);
	
	fit_df = pd.DataFrame(fits);		
	return(fit_df); 