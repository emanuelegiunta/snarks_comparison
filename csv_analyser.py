from csv_plotter import *

def ratio_vec(vec):
	if len(vec) < 2:
		raise CompatibilityError("Only one vector provided")

	for i in range(len(vec)):
		if len(vec[i]) != len(vec[0]):
			raise CompatibilityError("Vector of different length provided")

	out = []
	for j in range(1, len(vec)):
		# compute the component-wise ratio between vec[0] and vec[1]
		n = len(vec[0])
		out.append([round(vec[0][i]/vec[j][i], 2) for i in range(n)])

	return out

def test(filename):
	vec = csv_parse(filename)
	vec = vec[1:] # exclude the first line that is usually for the header
	rate = ratio_vec(vec)

	out = "For data in {}\n".format(filename)
	for rate_row in rate:
		r_min = min(rate_row)
		r_max = max(rate_row)
		r_last = rate_row[-1]

		rate_row.sort()
		r_med = rate_row[len(rate_row)/2]
		out += "  Ratio: min = {}, median = {}, max = {}, last = {}\n".format(r_min, r_med, r_max, r_last)

	print(out)

if __name__ == "__main__":
	argc = len(sys.argv)
	for i in range(1, argc):
		test(sys.argv[i])



