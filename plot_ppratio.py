import sys
import csv
import plotly.plotly as py      # Every function in this module will communicate with an external plotly server
from plotly.graph_objs import Scatter
from plotly.graph_objs import Data
from plotly.graph_objs import Layout
from plotly.graph_objs import Figure

#filename is 1st command line argument
csv_file = sys.argv[1]; 

rows = []

with open (csv_file, 'rU') as csvfile:
	wholefile = csv.reader(csvfile, delimiter=',', quotechar='|')
		
	for row in wholefile:
		
		rows.append([float(i) for i in row])
			

#even items are values		
#odd items are categories to create x axis
#converts every item to 
y_values = rows[0:][::2];

x_values = rows[1:][::2];

#marker styles
marker_y1 = dict(size = 12, symbol = 'square', color = 'rgb(102, 102, 102)', line=dict(color='rgb(0, 0, 0)', width=0.5))
marker_y2 = dict(size = 12, symbol = 'square', color = 'rgb(255, 255, 255)', line=dict(color='rgb(0, 0, 0)', width=0.5))
line = dict(dash = 'dot', width = '0.5', color='rgb(0, 0, 0)') ; 

#iteratively create traces
traces = []
for cell in range(len(y_values)):
	if x_values[cell][0] == 1:
		print('found 1')
		trace = Scatter(x=x_values[cell], y=y_values[cell], mode='lines+markers', marker=marker_y1);  
		traces.append(trace);
		
	if x_values[cell][0] == 2.5:
		print('found 3')
		trace = Scatter(x=x_values[cell], y=y_values[cell], mode='lines+markers', marker=marker_y2);  
		traces.append(trace);  

# (2) Make Data object 
data = Data(traces)  # (!) Data is list-like, must use [ ]

# (3) Make Layout object (Layout is dict-like)
# axis settings

y_axis_settings=dict(range = [0.8, 2.2], showline = True, title = "paired pulse ratio (40 msec)", titlefont=dict(family='Arial, sans-serif', size=20, color='rgb(0, 0, 0)'), tickfont=dict(family='Arial, sans-serif', size=18, color='rgb(0, 0, 0)'), showgrid = False, zeroline = False, linewidth = 0.5, color='rgb(0, 0, 0)') ; 
x_axis_settings=dict(range = [0, 5], showline = True, title = "", titlefont=dict(family='Arial, sans-serif', size=20, color='rgb(0, 0, 0)'), tickfont=dict(family='Arial, sans-serif', size=18, color='rgb(0, 0, 0)'), showgrid = False, zeroline = False, linewidth =0.5, color='rgb(0, 0, 0)') ; 

layout = Layout(title = '', yaxis = y_axis_settings, xaxis = x_axis_settings, width=475, height=575)

# (4) Make Figure object (Figure is dict-like)
fig = Figure(data=data, layout=layout) 

py.plot(fig, filename= csv_file.rstrip('.csv'), world_readable=False)







