import matplotlib.pyplot as plt
import math
import sys

class CompatibilityError(Exception):
	pass

def csv_parse_line(string):
	# remove the trailing \n
	string = string.replace("\n", "")
	
	# split with comma to get the values
	vec = string.split(",")
	vec = [ eval(x) for x in vec ]

	return vec

def vec_transpose(vec):
	assert len( {len(col) for col in vec} ) <= 1, "Transponing vector with different length lines"
	
	if vec == []:
		return vec
	else:
		n = len(vec[0])
		return [[line[i] for line in vec] for i in range(n)]

def csv_parse(file_name, transpose = True):
	f = open(file_name, "r")
	vec = []

	for line in f:
		vec.append(csv_parse_line(line))

	f.close()

	# we transpose so that lists of vec correspond to
	#  columns of file_name.csv unless transpose = False
	if transpose:
		vec = vec_transpose(vec)

	return vec

def _csv_plot(vec, ax):
	# debug
	assert len(vec) >= 2, "Not enough data to plot, expecting at least two columns"
	assert type(vec[0]) in (list, tuple), "Expected iterable of lists or tuples"
	assert len( {len(col) for col in vec} ) <= 1, "Plotting columns of different length, currently not supported"

	# get the x-axis values
	x = vec[0]
	# if the first entry is a string, use that for the x-label
	if type(x[0]) == str:
		ax.set_xlabel(x[0])
		x = x[1:]

	for i in range(1, len(vec)):
		# fancy coloring - just for plotting data form optimiser.py
		if i == 1:
			c = 'b'
		else:
			c = 'r'

		data = vec[i]
		legend_flag = False
		# if the first entry is a string, use that for the legend
		if type(data[0]) == str:
			data = data[1:]
			legend_flag = True

		line, = ax.plot(x, data, c = c)
		# if a legend was detected
		if legend_flag:
			line.set_label(vec[i][0])
			ax.legend()

def csv_plot(file_name, ax):
	assert file_name.split(".")[-1] == "csv", "Plotting non csv file"

	vec = csv_parse(file_name)
	_csv_plot(vec, ax)

	ax.set_ylim(ymin = 0)
	ax.set_xscale('log', basex = 2)
	ax.title.set_text(file_name)

if __name__ == "__main__":
	argc = len(sys.argv)

	n = argc - 1
	
	if n > 0:
		nrows = (n/1.6180)**(0.5)
		nrows = max(int(math.floor(nrows)), 1)

		ncols = (n/float(nrows))
		ncols = max(int(math.ceil(ncols)), 1)

		fig, axs = plt.subplots(ncols = ncols, nrows = nrows)
		
		for i, ax in enumerate(fig.axes):
			# arguments we need are from 1 to n
			csv_plot(sys.argv[i+1], ax)

		plt.show(fig)
