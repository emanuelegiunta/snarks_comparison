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
	vec = [ float(x) for x in vec ]

	return vec

def csv_parse(file_name):
	f = open(file_name, "r")
	vec = []

	for line in f:
		vec.append(csv_parse_line(line))

	f.close()
	return vec

def csv_plot_(vec, ax):
	if type(vec[0]) not in (list, tuple):
		raise CompatibilityError("csv_plot expects an iterable of lists or tuples")

	n = len(vec[0])
	for x in vec:
		if len(x) != n:
			raise CompatibilityError("csv_plot inputs don't have the same length")

	x = vec[0]
	for i in range(1, len(vec)):
		if i == 1:
			ax.plot(x, vec[i], c = 'r', lw = 2.0)
		else:
			ax.plot(x, vec[i], c = 'b', ls = "--")	

def csv_plot(file_name, ax):
	if file_name.split(".")[-1] != "csv":
		raise CompatibilityError("Attempted to plot a non-csv file")

	vec = csv_parse(file_name)
	csv_plot_(vec, ax)

	ax.set_ylim(ymin = 0)
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
