import stf





#TSeries08312017_1045_3339['voltage recording'].loc['Sweep0001'].loc[:,'time']

def plot_sweep(sweep_num, t_series_folder):
	
	str_sweep_num = str(sweep_num)
	
	y_vals = t_series_folder['voltage recording'].loc['Sweep'+str(str_sweep_num.zfill(4))].loc[:,'input 0']
	
	stf.new_window(y_vals) 
	
	time = t_series_folder['voltage recording'].loc['Sweep'+str(str_sweep_num.zfill(4))].loc[:,'input 0'].as_matrix()
	
	stf.set_sampling_interval((time[1]-time[0])/1000)
	
	return()
	
	