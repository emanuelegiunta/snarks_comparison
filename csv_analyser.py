from csv_plotter import *

def check_header(vec):
	return type(vec[0][0]) == str

def ratio_vec(vec):
	assert len(vec) >= 2, "Not enought data to make a comparison"
	assert len( {len(col) for col in vec} ) <= 1, "Columns have different length"

	# if there is a header in the csv
	#  ignore the first entry of each column
	if check_header(vec):
		m = 1
	else:
		m = 0

	out = []
	for j in range(1, len(vec)):
		# compute the component-wise ratio between vec[0] and vec[1]
		n = len(vec[0])
		out.append([round(vec[0][i]/float(vec[j][i]), 2) for i in range(m, n)])

	return out

def test(filename):
	vec = csv_parse(filename)
	vec = vec[1:] # exclude the first line that is usually for the number of constraints
	rate = ratio_vec(vec)

	out = "For data in {}\n".format(filename)
	for i, rate_row in enumerate(rate):
		r_min = min(rate_row)
		r_max = max(rate_row)
		r_last = rate_row[-1]

		if check_header(vec):
			row_name = vec[i+1][0]
		else:
			row_name = "column {:d}".format(i)

		rate_row.sort()
		r_med = rate_row[len(rate_row)//2]
		out += "  Ratio {:s}\t: min = {:.2f}, median = {:.2f}, max = {:.2f}, last = {:.2f}\n". \
			format(row_name, r_min, r_med, r_max, r_last)

	print(out)

if __name__ == "__main__":
	argc = len(sys.argv)
	for i in range(1, argc):
		test(sys.argv[i])



