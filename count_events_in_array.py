import sys
import numpy as np

def count_events_in_array(trace, threshold=0, up=True):
	""" Counts the number of events (e.g action potentials (AP)) in the current trace.
	Arguments:
	trace -- 1d numpy array
	dt -- sampling interval
	threshold -- (optional) detection threshold (default = 0).
	up -- (optional) True (default) will look for upward events, False downwards.

	Returns:
	An integer with the number of events.
	sample indicies of threshold crossings

	"""

	selection = trace
	# algorithm to detect events
	EventCounter,i = 0,0 # set counter and index to zero
	# list of sample points
	sample_points = []
	# choose comparator according to direction:
	if up:
		comp = lambda a, b: a > b
	else:
		comp = lambda a, b: a < b
	# run the loop
	while i<len(selection):

		if comp(float(selection[i]),float(threshold)):
			EventCounter +=1
			sample_points.append(i)

			i += 1 #skip set number of samples from threshold crossing
			
			while i<len(selection) and comp(selection[i],threshold):
				i+=1 # skip values if index in bounds AND until the value is below/above threshold again
		else:
			i+=1

	return (EventCounter, sample_points)