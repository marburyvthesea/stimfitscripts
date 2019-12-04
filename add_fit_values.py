
# coding: utf-8

# In[ ]:

import sys
add_to_path = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/wx-3.0-osx_cocoa', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks', '/Applications/MacPorts/stimfit.app/Contents/MacOS/../Frameworks/stimfit', '.', '', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Applications/MacPorts/stimfit.app/Contents/Frameworks/stimfit', '/Applications/MacPorts/stimfit.app/Contents/Frameworks', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/Library/Python/2.7/site-packages']
more_add_to_path = ['', '/Library/Python/2.7/site-packages/openpyxl-2.3.0-py2.7.egg', '/Library/Python/2.7/site-packages/et_xmlfile-1.0.1-py2.7.egg', '/Library/Python/2.7/site-packages/jdcal-1.2-py2.7.egg', '/Library/Python/2.7/site-packages/XlsxWriter-0.7.7-py2.7.egg', '/Library/Python/2.7/site-packages', '/Library/Python/2.7/site-packages/neurphys-0.0-py2.7.egg', '/Users/johnmarshall', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/analysis_with_excel', '/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/stimfitscripts', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC', '/Library/Python/2.7/site-packages'];
for path in add_to_path:
    sys.path.append(path);  
for path in more_add_to_path:
    sys.path.append(path);
import fitexpt
import stf
import pandas as pd
import numpy as np 
import matplotlib.pyplot as mplt


# In[ ]:

path_to_master = '/Users/johnmarshall/Documents/Analysis/RecordingData/masterabfspreadsheet.xlsx';
master_df = pd.read_excel(path_to_master);  


# In[ ]:

def get_file_paths(experiment_type_to_load):
    experiment_files = np.unique(master_df.loc[master_df['experiment type']==experiment_type_to_load]['associated experiment'].as_matrix());
    file_paths = [];
    for expt_file in experiment_files: 
        file_paths.append(master_df.loc[master_df['associated experiment']==expt_file]['directory'].as_matrix()[0]+'/'+expt_file);
    return(file_paths)

def get_raw_file_paths(experiment_type_to_load):
    experiment_files = np.array(master_df.loc[master_df['experiment type']==experiment_type_to_load]['file'].as_matrix());
    file_paths = [];
    for expt_file in experiment_files: 
        file_paths.append(master_df.loc[master_df['file']==expt_file]['directory'].as_matrix()[0]+'/'+expt_file);
    return(file_paths)   

def load_files_from_master(experiment_files_paths_series):
    files_dict = {}; 
    for fname in experiment_files_paths_series:  
    	csv_file = fname+'_stimfitanalysis.csv'; 
        print(csv_file); 
        try:
            files_dict[fname[-11:]] = pd.read_csv(csv_file); 
            #replace loaded indicies with time values
            indicies_length = files_dict[fname[-11:]]['Unnamed: 0'].shape[0];
            time = np.arange(0,float(indicies_length)*(float(1)/float(3)), float(1)/float(3));
            files_dict[fname[-11:]]['Unnamed: 0'] = time ; 
        except IOError:
            print('could not find', csv_file)
            pass
    
    files_panel = pd.Panel(files_dict);
    mean_panel = files_panel.mean(axis=0);
    return(files_panel, mean_panel)


# For fitting capacitance transients from .abf files
def batch_fit(file_list, fit_search_params, pulse_length, fit_function):

	output_list = []
	for fname in file_list:
		print(fname);
		stf.file_open(str(fname)); 
		fit_df_exp = fitexpt.fit_experiment(fit_search_params, pulse_length, fit_function);
		
		fit_df_exp.to_csv(fname+'_fit_results.csv', sep=',');
		output_list.append(fname+'_fit_results.csv'); 
	
	
	
	
	return(output_list)
# In[ ]:

def write_to_excel_with_charts_formulas(directory, file_name, panel):
    writer = pd.ExcelWriter(directory+file_name, engine='xlsxwriter'); 
    for df in panel:
        frame_to_write = panel[df];
        frame_to_write.to_excel(writer, sheet_name=df); 
        workbook = writer.book;
        worksheet = writer.sheets[df];
        charts_to_write = ['EPSC1', 'EPSC2', 'Rs(Mohm)', 'Ri', 'Capacitance(pA)'];
        columns_with_data = {'EPSC1':7,'EPSC2':8, 'Rs(Mohm)':3, 'Ri':5, 'Capacitance(pA)':6}
        position_ = [4, 20, 36, 52, 68];
        for chart_name, row_position in zip(charts_to_write, position_): 
            chart = workbook.add_chart({'type': 'line'}); 
            chart.add_series({'values':[df, 1,columns_with_data[chart_name],125,columns_with_data[chart_name]],
                              'categories': [df, 1, 1, 125, 1]});
            chart.set_x_axis({'name': 'time(min)', 'position_axis': 'on_tick', 'interval_unit': 30, 'interval_tick': 30}); 
            chart.set_y_axis({'name': chart_name, 'major_gridlines': {'visible': False}}); 
            # Turn off chart legend. It is on by default in Excel.
            chart.set_legend({'position': 'none'}); 
            # Insert the chart into the worksheet.
            worksheet.insert_chart('R' + str(row_position), chart);
    #write formula headers
        bold = workbook.add_format({'bold': 1})
        columns = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2'];
        header_positions = ['K1','L1','M1','N1','O1','P1','Q1'];
        for column, position in zip(columns, header_positions):
            worksheet.write(position, column, bold);
    #write formulas for baseline means
        worksheet.write('J2', 'Baseline Average', bold);
        columns_to_average = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2'];
        data_ranges = {'CapTransPeak(pA)':'C2:C32','Rs(Mohm)':'D2:D32','SS(pA)':'E2:E32','Ri':'F2:F32','Capacitance(pF)':'G2:G32','EPSC1':'H2:H32','EPSC2':'I2:I32'};
        formula_position = ['K2', 'L2', 'M2', 'N2', 'O2', 'P2', 'Q2']
        for chart_name, position in zip(columns_to_average, formula_position):
            worksheet.write_formula(str(position), '{=AVERAGE('+str(data_ranges[chart_name])+')}'); 
    #write formulas for KA application means
        worksheet.write('J4', '10-15 min KA Average', bold);
        columns_to_average = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2'];
        data_ranges = {'CapTransPeak(pA)':'C62:C77','Rs(Mohm)':'D62:D77','SS(pA)':'E62:E77','Ri':'F62:F77','Capacitance(pF)':'G62:G77','EPSC1':'H62:H77','EPSC2':'I62:I77'};
        formula_position = ['K4', 'L4', 'M4', 'N4', 'O4', 'P4', 'Q4']
        for chart_name, position in zip(columns_to_average, formula_position):
            worksheet.write_formula(str(position), '{=AVERAGE('+str(data_ranges[chart_name])+')}'); 
    #write formulas for wash means
        worksheet.write('J6', '10-15 min Wash Average', bold);
        columns_to_average = ['CapTransPeak(pA)','Rs(Mohm)','SS(pA)','Ri','Capacitance(pF)','EPSC1','EPSC2'];
        data_ranges = {'CapTransPeak(pA)':'C92:C122','Rs(Mohm)':'D92:D122','SS(pA)':'E92:E122','Ri':'F92:F122','Capacitance(pF)':'G92:G122','EPSC1':'H92:H122','EPSC2':'I92:I122'};
        formula_position = ['K6', 'L6', 'M6', 'N6', 'O6', 'P6', 'Q6']
        for chart_name, position in zip(columns_to_average, formula_position):
            worksheet.write_formula(str(position), '{=AVERAGE('+str(data_ranges[chart_name])+')}'); 
    
    writer.save()
    return()



