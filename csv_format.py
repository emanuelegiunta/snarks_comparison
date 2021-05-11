from csv_plotter import *

def _csv_row_to_col(vec):
	wec = []
	for i in range(len(vec[0])):
		wec.append([row[i] for row in vec])

	return wec

def _discard_data(vec):
	return vec[2:]

def _bit_to_kb(vec):
	# rescale only the data column
	for i in range(1, len(vec[0])):
		for row in vec:
			row[i] = row[i]/(8 * 10**6)

def _log2_to_linear(vec):
	# rescale only the first column
	for row in vec:
		row[0] = 2**row[0]

def vec_to_csv(file, vec):
	#werite a vector into an open csv file
	file.write(",".join([str(a) for a in vec]) + "\n")

def csv_row_to_col(file_name):
	print("Adjusting {:s}".format(file_name))
	vec = csv_parse(file_name)
	wec = _csv_row_to_col(vec)
	wec = _discard_data(wec)

	_bit_to_kb(wec)
	_log2_to_linear(wec)

	header = [["x"] + ["data{:d}".format(i) for i in range(1, len(vec))]]
	wec = header + wec

	f = open("New_"+file_name, "w")
	for row in wec:
		vec_to_csv(f, row)
	f.close()


if __name__ == "__main__":
	argc = len(sys.argv)
	for i in range(1, argc):
		csv_row_to_col(sys.argv[i])