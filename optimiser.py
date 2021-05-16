import numpy as np
import numpy.random as rnd
import math
import matplotlib.pyplot as plt
import sys 

# - - Constants - - #
# constants used by BSC
MKT_ITER = 1				# recommended 2000	
MKT_ESTIMATE_ITER = 20		# recommended 20
MTK_FAST_FLAG = True 		# recommended True

# constants used by ligero.
LGR_MAX_SAMPLE = 150		# recommended 150
LGR_MAX_DOMAIN_DIM = 6		# recommended 6

# constants used by FRI
FRI_ITER = 1					# recommended 1000
FRI_MAX_DOMAIN_DIM = 6			# recommended 6
FRI_MAX_LOCALIZATION_NUM = 4	# recommended 4

class CompatibilityError(Exception): #de
	pass

class InputError(Exception):
	pass

# - - - - - - - - - - - - - - - - - - - #
#	basic math and generic functions	#
# - - - - - - - - - - - - - - - - - - - #

def comb(n,m):
	#binomial coefficient
	out = 1
	for i in range(m):
		out *= (n - i)/float(m - i)
	return out

def combr(n, m, k):
	#ratio of comb(n,k)/comb(m,k)
	out = 1
	for i in range(k):
		out *= (n - i)/float(m - i)
	return out

def print_cost(name, bit_cost):
	if VERBOSE_FLAG:
		print((name + ":\t%10d bit%10d B%10d KB") % (bit_cost, bit_cost/8, bit_cost/8000) )

def print_verbose_cost(name, bit_cost):
	if VERY_VERBOSE_FLAG:
		print_cost(name, bit_cost)

def print_message(string):
	if VERBOSE_FLAG:
		print(string)

def print_verbose_message(string):
	if VERY_VERBOSE_FLAG and VERBOSE_FLAG:
		print(string)

def print_separator(length = 36, indentation_length = 2, very_verbose = False):
	# print a line to separate data
	#  length 				: number of caracter printed
	#  indentation_length 	: number of spaces before the line
	#  very_verbose 		: if True check the VERY_VERBOSE_FLAG
	if (not very_verbose) or VERY_VERBOSE_FLAG:
		print("\n" + " "*indentation_length + "="*length + "\n")

def vec_to_csv(file, vec):
	#werite a vector into an open csv file
	file.write(",".join([str(a) for a in vec]) + "\n")

def str_bold(s):
	# print bold text, need the terminal to support ANSI
	return "\033[1m{:s}\033[0m".format(s)

def str_underline(s):
	# print underlined text, need the terminal to support ANSI
	return "\033[4m{:s}\033[0m".format(s)

def str_synopsis():
	out  = ""
	out += str_bold("SYNOPSIS: ") + "{}\n".format(sys.argv[0])
	out += "\t[-l] [-tl] [-a] [-ta] [-sp {arg}] [-ld {arg}] [-hd {arg}]\n"
	out += "\t[-h] [-p] [-r {args}] [-fd {arg}] [-f {arg}]\n"
	out += "\t[-D {arg}] [-v] [-vv] [-help]\n"
	out += "\n"
	
	out = out.format(arg = str_underline("1 argument"), args = str_underline("2 arguments"))
	return out

def str_help():
	out  = ""
	out += str_bold("Parameter optimiser for Aurora [EC:BCRSVW19] and Ligero [CCS:AHIV17]\n\n")
	out += str_synopsis()

	out += str_bold("COMMAND LINE DESCRIPTION:\n")
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-l, --ligero"))
	out += " Runs a comparison between standard and optimised ligero\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-a, --aurora"))
	out += " Runs a comparison between standard and optimised aurora\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-sp, --security_parameter"), str_underline("(1 argument)"))
	out += " Set the security parameter\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-ld, --lowest_dimension"), str_underline("(1 argument)"))
	out += " Set the minimum size of tested constraints, inclusive\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-hd, --highest_dimension"), str_underline("(1 argument)"))
	out += " Set the maximum size of tested constrants, inclusive\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-h, --heuristic"))
	out += " Uses optimistic soundness bounds\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-p, --proven"))
	out += " Uses proven soundness bounds\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-r, --rmfe"), str_underline("(2 arguments)"))
	out += " Set the reverse multiplication-friendly embedding (k,m) used by the optimised version\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-fd, --field_dimension"), str_underline("(1 argument)"))
	out += " Fix the field dimension over F_2 for the non optimised version\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-f, --file_name"), str_underline("(1 argument)"))
	out += " Set the file name. If not given derive the file name from the parameters\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-D, --directory"), str_underline("(1 argument)"))
	out += " Set the directory. as no check is performed on the given path can also be used to attach a prefix to the file name\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-v, --verbose"))
	out += " Print some information on screen\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-vv, --very_verbose"))
	out += " Print more informations on screen\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-help"))
	out += " Prints the following message\n"
	out += "\n"

	return out

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#	monte carlo cost evalution for merkel-trees opening   #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Faster Version

def mkt_estimate_cost(leaf_num, path_num, fast = MTK_FAST_FLAG):
	if fast:
		out = 0
		h = int(math.ceil(np.log2(leaf_num))) #estimate the height
		for i in range(1, h+1):
			width = 2.0**i
			#add the expected number of hash. Assumes collisions are possibile
			#for an explanation of this formula see libiop (which file?)
			out += path_num * ((width - 1)/width)**(path_num - 1) 

		out += path_num
		return int(round(out))
	else: #return a monte-carlo evaluation over less iterations
		return mkt_avg_cost(leaf_num, path_num, iter_num = MKT_ESTIMATE_ITER)

def mkt_single_cost(leaf_num, path_num):
	# assuming no collision is allowed, we sample a set of size path_num
	leaves_set = set() 	# Empty set

	#create a random set of path_num elements in range(leaf_num)
	while len(leaves_set) != path_num:
			delta = path_num - len(leaves_set)
			sample = set(rnd.randint(0, leaf_num, delta))
			leaves_set = leaves_set.union(sample)

	out = len(leaves_set)
	size_new = out

	k = int(math.ceil(np.log2(leaf_num)))
	f = lambda i, n : n % (1 << (k - i))
	i = 0
	
	# Count the number of hashes required per layer
	for i in range(1, k + 1):
		size_old = size_new
		leaves_set = {f(i,n) for n in leaves_set}
		size_new = len(leaves_set)

		# The number of node we need to count is the number of those already
		# in the previous layer (precisely their borthers) minus the number
		# of node whose brother was already present in the previous layer
		out += size_old - 2*(size_old - size_new)

	return out

def mkt_avg_cost(leaf_num, path_num, iter_num = MKT_ITER):
	out = 0.0

	for i in range(iter_num):
		out += mkt_single_cost(leaf_num, path_num)

	return out / iter_num

# - - - - - - - - - - - - - - - - - - - #
#	BCS transform communication costs	#
# - - - - - - - - - - - - - - - - - - - #

def com_BCS(hash_size, vec_alphabet_size, rounds, vec_oracle_length, vec_query, estimate = False):
	#hash_size : bit lambdaength of the hash output used (128/256)
	#vec_alphabet_size : bit size of the alphabet used for the oracles (i.e. size of the field)
	#rounds : number or rounds
	#vec_oracle_length : length of the oracles sent
	#vec_query : number of queries to the i-th oracle

	#adapt single values for vec_oracle_length and vec_query:
	if type(vec_oracle_length) != list:
		vec_oracle_length = [vec_oracle_length]

	if type(vec_query) != list:
		vec_query = [vec_query]

	if type(vec_alphabet_size) != list:
		vec_alphabet_size = [vec_alphabet_size]

	assert len(vec_query) == len(vec_oracle_length), \
		"BSC, length of queries ({:d}) and oracles ({:d}) differ".format(len(vec_query), len(vec_oracle_length))

	assert len(vec_query) == len(vec_alphabet_size), \
		"BSC, length of queries ({:d}) and alphabet ({:d}) differ".format(len(vec_query), len(vec_alphabet_size))

	if not estimate:
		#compute the number of mkt hashes required from Monte-Carlo evaluation
		mkt_opening = hash_size*sum(mkt_avg_cost(vec_oracle_length[i], vec_query[i]) for i in range(len(vec_query)))
		mkt_opening = int(math.ceil(mkt_opening))
	else:
		#compute the number of mkt hashes required from an estimate of the expected value, faster
		mkt_opening = hash_size*sum(mkt_estimate_cost(vec_oracle_length[i], vec_query[i]) for i in range(len(vec_query)))

	mkt_opened_elements = sum(vec_alphabet_size[i] * vec_query[i] for i in range(len(vec_query))) #cost of the opened items
	
	mkt_round_hash = hash_size*(rounds + 1) #cost of the seed for verifiers replies

	#print_cost("hashed path",mkt_opening)
	#print_cost("oracle answer", mkt_opened_elements)
	#print_cost("round hashes", mkt_round_hash)

	out = mkt_opening + mkt_opened_elements + hash_size*(rounds + 1) 
	return out

# - - - - - - - - - - - - - - - #
#	Ligero for binary circuits	#
# - - - - - - - - - - - - - - - #

class LIGERO_parameters:
	def __init__(self, n, m, fd, sp, snd_type, protocol_type, rmfe = None):

		self.protocol_type = protocol_type	# "standard" or "optimised"

		self.variables = n 					# variables
		self.constraints = m 				# constraints
		self.field_dim = fd 				# log2(|F|)	

		if self.protocol_type == "optimised":
			# debug - the rmfe parameter has to be given
			assert rmfe != None, "LIGERO, rmfe not given in optimised version"
			# - - - - - end of debug - - - - - #

			self.field_dim = rmfe[0]
			self.rmfe = rmfe[1]				# Warning : this is ambiguous

		# checking if the size of the field is fixed. A not given
		#  value implies that the field size can be adjusted to minimize
		#  communication costs.
		if self.field_dim == None:
			self.field_type = "variable"
			self.field_size = None						
		else:
			self.field_type = "fixed"
			self.field_size = 2**self.field_dim

		# setting the number of variables/constraints
		if self.protocol_type == "standard":
			# for each variable we need the constraint x^2 = x
			self.constraints += self.variables

		elif self.protocol_type == "optimised":
			#variables and constrains scale down of a rmfe factor
			self.variables = int(math.ceil(self.variables/float(self.rmfe)))
			self.constraints = int(math.ceil(self.constraints/float(self.rmfe)))

		self.hash_size = 2*sp
		self.security_parameter = sp 	 	# kappa

		if self.protocol_type == "standard":
			self.interactive_soundness_error = sp + 1		# one term
		elif self.protocol_type == "optimised":
			self.interactive_soundness_error = sp + 3		# four terms

		self.query_soundness_error = sp + 1					# one term

		self.domain_dim = None				# n
		self.domain_size = None 			# 2^n 		**
		self.degree = None					# k 		
		self.proximity_parameter = None 	# e 		**
		self.queries = None 				# t 		**
		self.factor_l = None				# l 		
		self.factor_m1 = None				# m1 		** number of codewords to fill w
		self.factor_m2 = None				# m2 		** number of codewords to fill x, y, z, [t]
		self.interactive_repetitions = None	# sigma 	**

		self.snd_type = snd_type

		# constants
		if self.protocol_type == "standard":
			self.LGR_MIN_L = 1/10.0
			self.LGR_MAX_L = 2/3.0

		elif self.protocol_type == "optimised":
			self.LGR_MIN_L = 1/6.0
			self.LGR_MAX_L = 2.5

	def check_empty_entries(self, var = [], exclude = False):
		# check that all the variables beside those whose name
		#  is in the list "all_but" are not set to None
		# 
		# It returns the name of the empty variable. If no empty
		#  variable is found, it returns None

		for key, value in vars(self).items():
			if value == None and (key in var or exclude) and not (key in var and exclude):
				return key
		else:
			return None

	def __str__(self):
		out = "\n- - - - LGR pramaters - - - -\n"

		out += "protocol type:\t\t%s\n" % self.protocol_type
		out += "soundness type:\t\t%s\n" % str(self.snd_type)

		out += "variables:\t\t%d\n" % self.variables
		out += "constraints:\t\t%d\n" % self.constraints
		
		out += "field dimension:\t%s\n" % str(self.field_dim)
		out += "field type:\t\t%s\n" % str(self.field_type)

		if self.protocol_type == "optimised":
			out += "RMFE parameters:\t({},{})\n".format(self.rmfe, self.field_dim)

		out += "query soundness:\t%d\n" % self.query_soundness_error
		out += "interactive soundnes:\t%d\n" % self.interactive_soundness_error

		out += "RS block length:\t%s\n" % str(self.domain_dim)
		out += "RS degree:\t\t%s\n" % str(self.degree)
		out += "proximity parameter:\t%s\n" % str(self.proximity_parameter)
		out += "queries:\t\t%s\n" % str(self.queries)
		out += "parameter l:\t\t%s\n" % str(self.factor_l)
		out += "parameter m1:\t\t%s\n" % str(self.factor_m1)
		out += "parameter m2:\t\t%s\n" % str(self.factor_m2)
		out += "retepetitions:\t\t%s\n" % str(self.interactive_repetitions)

		out += "- - - - - - - - - - - - - - -"
		return out

	def find_queries(self):
		# Returns the number of queries required to have a query soundness error smaller
		#  that 2^{-qsp}

		# debug - check that all the required variables are non empty
		assert self.check_empty_entries(var = ("domain_size", "proximity_parameter", "degree", \
			"query_soundness_error")) == None, "LIGERO, empty required variable"
		# - - - - - - - - - - debug - - - - - - - - - - #

		# we find t through a binary search
		t_min = 0
		t_max = 1

		n = self.domain_size
		e = self.proximity_parameter
		k = self.degree
		err = 2.0**(-self.query_soundness_error)

		# f(t) returns the query soundness error with t queries
		f = lambda t : max(combr(n - e - 1, n, t), combr(e + 2*k - 2, n, t))

		# we rule out the case in which soundness cannot be obtained even
		#  if all points of the inputs are queried
		if f(n) >= err:
			return float('inf')

		#first we look for an upper bound
		while(f(t_max) >= err):
			t_max = 2*t_max
		
		t_min = t_max/2
		t_max = t_max

		# next we run the binary search to minimize f
		while(t_max - t_min > 1):
			if f((t_min + t_max)/2) >= err:
				t_min = (t_min + t_max)/2

			else:
				t_max = (t_min + t_max)/2

		return t_max

	def find_field_dim(self):
		# IDEA - significant efficiency improvement (expected 1.5x) could come
		#  using 3 different repetition parameters sigma_LDT, sigma_LIN, sigma_ROW.
		#  (TODO)
		
		# debug - check that all the required variables are not-empty
		assert self.check_empty_entries(var = ("proximity_parameter", \
			"interactive_soundness_error", "snd_type", "domain_dim")) == None, \
			"LIGERO, empty required variable".format()
		# - - - - - - - - - - end of debug - - - - - - - - - - #

		pp = np.log2(self.proximity_parameter + 1)
		isp = self.interactive_soundness_error

		if self.snd_type == "proven":
			# assuming an interactive soundness error of ((e + 1)/|F|)^sigma
			sigma = lambda f : int(math.ceil( isp/(f - pp) ))
		
		elif self.snd_type == "heuristic":
			# assuming an interactive soundness error of (e + 1)/|F|^sigma
			sigma = lambda f : int(math.ceil( (isp + pp)/f ))		
		
		# cost(f) returns the cost of executing LIGERO with the current
		#  parameters excluding the mkt in BSC - it therefore include
		#  the costs from the queries (the replies without the accepting
		#  path)
		cost = lambda f : self.find_other_cost(sigma = sigma(f), f = f, query_flag = True)

		cost_min = float('inf')
		f_min = None
		cost_tmp = None

		# we search element by element the best f scanning
		#  the entire range. Remark that isp + pp + 1 is the value that makes
		#  the number of repetitions 1 in the worst calse - further 
		#  increasing f does not provide improvements
		for f in range(self.domain_dim, int(math.ceil(isp + pp + 1))):

			cost_tmp = cost(f)
			if cost_tmp < cost_min:
				cost_min = cost_tmp
				f_min = f

		return f_min

	def find_other_cost(self, sigma = None, f = None, query_flag = False):
		# Returns the cost of sending plain messages. if query_flag
		#  is True it also considers the cost of opening oracles
		#  without the accepting path
		#
		# sigma and f are provided to make "speculative" calls to 
		#  this function without actually modifying self.field_dim
		#  or self.interactive_repetitions

		# debug - check for empty required entries
		assert self.check_empty_entries(var = ("domain_size", "degree", "factor_l", "protocol_tpye", \
			"queries")) == None, "LIGERO, empty required entries"
		# - - - - - - - - - - end of debug - - - - - - - - - - #

		if sigma == None:
			sigma = self.interactive_repetitions
		if f == None:
			f = self.field_dim

		n = self.domain_size
		k = self.degree
		l = self.factor_l
		#sigma = self.interactive_repetitions
		#f = self.field_dim

		# we set 3*(k + l - 1) as we only perform 3 lincheck in the R1CS
		#  notice that we perform the same number of lincheck/rowchecks
		out = f * sigma * (n + 3*(k + l - 1) + (2*k - 1))

		if self.protocol_type == "optimised":
			# in the optimised version we send two vectors of len isp
			out += 2 * f * self.interactive_soundness_error

		if query_flag:
			out += self.queries * self.find_alphabet_size(sigma = sigma, f = f)

		return(out)

	def find_rounds(self):
		if self.protocol_type == "standard":
			return 2
		elif self.protocol_type == "optimised":
			return 3

	def find_alphabet_size(self, sigma = None, f = None):
		# return the alphabet size of the (column compressed) oracles
		#  used in BSC.
		# 
		# sigma and f are provided to make "speculative" calls. if not
		#  provided we set them as the associated values in self

		# debug - check for empty entries
		assert self.check_empty_entries(var = ("protocol_type", "factor_m1", "factor_m2")) == None, \
			"LIGERO, empty required variable"

		if sigma == None:
			# debug - check for empty entries
			assert self.check_empty_entries(var = ("interactive_repetitions")) == None, \
				"LIGERO, empty required variable"
			# - - - - - end of debug - - - - - #
			sigma = self.interactive_repetitions

		if f == None:
			# debug - check for empty entries
			assert self.check_empty_entries(var = ("field_dim")) == None, \
				"LIGERO, empty required variable"
			# - - - - end of debug - - - - - #
			f = self.field_dim

		if self.protocol_type == "standard":
			# oracles: w (m1), x, y, z (m2)
			# masking terms: sigma for 1 LDT, 3 Lin, 1 Row
			out = (self.factor_m1 + 3*self.factor_m2 + 5 * sigma) * f
		elif self.protocol_type == "optimised":
			# oracles: w (m1), x, y, z, t (m2)
			# masking terms: sigma for 1 LDT, 3 Lin_h, 1 Row
			out = (self.factor_m1 + 4*self.factor_m2 + 5 * sigma) * f

			# debug - check for empty entries
			assert self.check_empty_entries(var = ("interactive_soundness_error", "factor_l")) == None, \
				"LIGERO, empty required variable"
			# - - - - - end of debug - - - - - #

			# masking terms for the modular lincheck
			out += int(math.ceil((3 * self.interactive_soundness_error)/self.factor_l)) * f

		return out

	def find_l_0(self):
		# Returns the base point l_0 s.t. the protocol looks for
		#  the best cost in [LGR_MIN_L * l_0, LGR_MAX_L * l_0]
		
		# debug - check for empty entries
		assert self.check_empty_entries(var = ("variables", "constraints", "security_parameter")) == None, \
			"LIGERO, empty required variable"
		# - - - - - end of debug - - - - - #

		l_0 = (max(self.variables, self.constraints)*float(self.security_parameter))**(0.5)

		return(l_0)

	def set_factor_m(self):
		# m1 = ceil(n/l) so that l m1 >= n
		self.factor_m1 = self.variables / float(self.factor_l)
		self.factor_m1 = int(math.ceil(self.factor_m1))

		# m2 = ceil(m/l) so that l m2 >= n
		self.factor_m2 = self.constraints / float(self.factor_l)
		self.factor_m2 = int(math.ceil(self.factor_m2))

	def complete(self):
		# Using degree, domain_dim and factor_l set the other parameters
		#  according to [AHIV17]

		# debug - check for empty variables
		assert self.check_empty_entries(var = ("domain_dim", "factor_l", "degree")) == None, \
			"LIGERO empty required variable"
		# - - - - - - - - - - end of debug - - - - - - - - - - #

		# number of points in the RS domain
		self.domain_size = 2**self.domain_dim

		# e = floor((n - 2k + 1)/4)
		# if we assume strong_ligero 4 -> 3
		self.proximity_parameter = (self.domain_size - 2*self.degree + 1)/3.0
		self.proximity_parameter = int(math.floor(self.proximity_parameter))

		# set the factors m1 and m2
		self.set_factor_m()

		# t is the minimum value that makes the query soundness error smaller
		#  than 2^(-qsp)
		self.queries = self.find_queries()
		
		if self.queries != float('inf') and self.field_type == "variable":
			# if the parameters allow a small query soundness and the field
			#  size is not fixed we find the field size that minimises costs
			self.field_dim = self.find_field_dim()
			self.field_size = 2**self.field_dim
		
		elif self.field_type == "variable":
			# if the parameters don't allow a small query soundness error
			#  we could exit the function. Instead we choose to set field
			#  dimension and other parameters to a recognizable state for debug
			self.field_dim = 1
			self.field_size = 1

		if self.queries != float('inf'):
			# if the field does allow for a small query soundness error

			if self.snd_type == "proven":
				# sigma is the minimum value that makes the interactive soundness error
				#  smaller than 2^(-isp)
				self.interactive_repetitions = \
					(self.interactive_soundness_error)/(self.field_dim - np.log2(self.proximity_parameter + 1))
				self.interactive_repetitions = int(math.ceil(self.interactive_repetitions))
			
			elif self.snd_type == "heuristic":
				self.interactive_repetitions = \
					(self.interactive_soundness_error + np.log2(self.proximity_parameter + 1))/float(self.field_dim)
				self.interactive_repetitions = int(math.ceil(self.interactive_repetitions))

		# debug - check that computation was correct
		assert self.queries >= 1, "LIGERO, non positive queries"
		# - - - - - end of debug - - - - - #

	def estimate_cost(self, estimate = True):
		# return an estimate of Ligero's cost
		#  with the current parameters. if estimate
		#  is True BSC will be faster and less precise

		# debug - check for empty entries
		assert self.check_empty_entries(var = ("domain_size", "queries", "hash_size")) == None, \
			"LIGERO, empty required variable"
		# - - - - - - - - - end of debug - - - - - - - - - - #

		rounds = self.find_rounds()
		alphabet_size = self.find_alphabet_size()
		oracle_length = self.domain_size
		queries = self.queries

		bsc_cost = com_BCS(self.hash_size, alphabet_size, rounds, oracle_length, queries, estimate)

		other_cost = self.find_other_cost()

		print_verbose_cost("\tBSC cost", bsc_cost)
		print_verbose_cost("\tOther cost", other_cost)

		return bsc_cost + other_cost

	def avg_cost(self):
		# return the average cost of running ligero with 
		#  current parameters - slow but precise
		return self.estimate_cost(estimate = False)

	def optimize(self, plot = False):
		# return self with the best parameters based on the witness size, and the
		#  security parameter. If plot is True return a plot of the cost as a
		#  function of l
		#
		# the only three free parameters are k, n, m with the constraints
		#  k >= l + t
		#  dim_n > log2(k)

		l_0 = self.find_l_0()
		l_min = int(math.ceil(l_0 * self.LGR_MIN_L))
		l_max = int(math.ceil(l_0 * self.LGR_MAX_L))
		l_step = int(math.ceil((l_max - l_min)/float(LGR_MAX_SAMPLE)))

		final_cost = float('inf')			# resulting cost
		final_tuple = None					# tuple associated with the best cost

		if plot:							# if plotting cost as a function of l
			vec_l = []						#  init a list for the checked l-values
			vec_cost = []					#  init a list for the associated cost

		# we make m vary in a constant fixed range.
		for self.factor_l in range(l_min, l_max, l_step):

			if plot:
				local_cost = float('inf')

			# k = l is the minimum value for the 
			#  degree of the RS code. In general k >= l + t
			self.degree = self.factor_l 	

			# list of all the domain values to check. 
			#  the base value is log2(2k + 1) to make 
			#  n >= 2k + 1.
			domain_values = [int(math.ceil(np.log2(2*self.degree + 1))) + i for i in range(LGR_MAX_DOMAIN_DIM)]
			domain_values_tmp = []

			while(domain_values):								# while there are still values to check
				min_degree_bound = float("inf")

				for self.domain_dim in domain_values:

					self.complete()								# set the other parameters as functions of m,k,n
					if self.degree - self.factor_l >= self.queries:
						print_verbose_message("\tTesting (l, k, log2n, delta) = ({:d}, {:d}, {:2d}, {})".format(
							self.factor_l, self.degree, self.domain_dim, self.degree - self.factor_l - self.queries))

						cost_tmp = self.estimate_cost()			# estimate costing with curretn m,k,n
						domain_values.remove(self.domain_dim)	# removing the current domain dim

						print_verbose_cost("\tProof size", cost_tmp)
						print_verbose_message("")

						if cost_tmp < final_cost:				# if the current cost is better than the current minimum
							final_cost = cost_tmp				# record it
							final_tuple = (self.factor_l, self.degree, self.domain_dim)

						if plot and cost_tmp < local_cost:		# if I'm plotting costs as a function of l
							local_cost = cost_tmp				#  i store the local minimum wrt the current l

					elif self.queries != float('inf'):
						domain_values_tmp.append(self.domain_dim)
						min_degree_bound = min(min_degree_bound, self.factor_l + self.queries)

				domain_values, domain_values_tmp = domain_values_tmp, []	
				self.degree = min_degree_bound

			# after all possibilities for the current l are tested
			#  if I'm plotting, i store the values
			if plot:
				vec_l.append(self.factor_l)
				vec_cost.append(local_cost)

		# after scanning for all n, k, m:
		self.factor_l, self.degree, self.domain_dim = final_tuple	# get the best parameters
		self.complete()												#  adjust other parameters

		if plot:
			plt.plot(vec_l, vec_cost)
			plt.show()

	def optimal_cost(self, plot = False):
		# optimise all parameters and return the total proof size
		self.optimize(plot = plot)
		return self.avg_cost()

def LIGERO_test(dim = 18, sp = 128, fd = None, rmfe = None, snd_type = "proven", protocol_type = "standard", plot = False):
	print_message("testing {:s} ligero, n = {:2d}".format(protocol_type, dim))
	ligero_p = LIGERO_parameters(2**dim, 2**dim, sp, fd = fd, rmfe = rmfe,\
		snd_type = snd_type, protocol_type = protocol_type)
	cost = ligero_p.optimal_cost(plot = plot)

	print_message(str(ligero_p))
	return cost

# - - - - - - - - - - - - - - - - - #
#	FRI soundness error and costs	#
# - - - - - - - - - - - - - - - - - #

class FRI_parameters:
	def __init__(self, n, m, fd, isp, qsp, h, snd_type = "proven", other_oracles = []):

		self.variables = n
		self.constraints = m
		self.field_dim = fd
		self.field_size = 2**fd #redundant
		self.hash_size = h

		self.query_soundness_error = qsp		
		self.interactive_soundness_error = isp	# Note: currently the isp is ignored

		self.domain_dim = None
		self.domain_size = None	#redundant
		self.rate = None
		self.query_bound = None
		self.query_repetition = None 
		self.localization_number = None
		self.relative_distance = None

		# List of other oracles containing the number of codewords per round i.e.
		# [n_1, ..., n_r] with n_i an integer. This is converted in actual oracle data
		# i.e. alphabet size, proof length and queries by the methods
		# 	coset_hash_alphabet()
		#	coset_hash_length()
		#	coset_hash_queries()
		self.other_oracles = other_oracles

		# Either "proven" or "heuristic"
		self.snd_type = snd_type

	def check_empty_entries(self):
		if self.domain_dim == None:
			return False

		if self.domain_size == None:
			return False

		if self.rate == None:
			return False

		if self.query_bound == None:
			return False

		if self.query_repetition == None:
			return False

		if self.localization_number == None:
			return False

		if self.relative_distance == None:
			return False

		return True

	def __str__(self):
		out = "\n- - - - FRI pramaters - - - -\n"

		out += "variables:\t\t%d\n" % self.variables
		out += "constraints:\t\t%d\n" % self.constraints
		out += "field dimension:\t%d\n" % self.field_dim
		out += "hash output size:\t%d\n" % self.hash_size
		
		out += "query soundness:\t%d\n" % self.query_soundness_error
		out += "interactive soundness:\t%d\n" % self.interactive_soundness_error

		out += "RS block lenght:\t%s\n" % str(self.domain_dim)
		out += "RS code rate:\t\t%s\n" % str(self.rate)
		out += "queries bound:\t\t%s\n" % str(self.query_bound)
		out += "queries repetitions:\t%s\n" % str(self.query_repetition)
		out += "localization number:\t%s\n" % str(self.localization_number)
		out += "relative distance:\t%s\n" % str(self.relative_distance)
		out += "soundness type:\t\t%s\n" % str(self.snd_type)

		out += "- - - - - - - - - - - - - - -"

		return out

	def estimate_queries(self, N):
		# Montecarlo simulation to estimate the number of
		# distinct queries to coset-hashed oracles in layer L^(i)
		# with i such that |L^(i)| = N
		estimated_queries = 0.0

		for j in range(FRI_ITER):									
			sample = rnd.randint(N, size = self.query_repetition)	#sample l points
			estimated_queries += len(set(sample))					#consider only the distinct ones

		estimated_queries = int(math.ceil(estimated_queries/float(FRI_ITER)))
		return estimated_queries

	def coset_hash_alphabet(self):
		# Set the alphabet size for the other oracles (leaves are grouped in sets of 2^eta)
		return [int(math.ceil( n * self.field_dim * (2.0**self.localization_number) )) for n in self.other_oracles]

	def coset_hash_length(self):
		# Set the proof length for the other oracles (leaves are grouped in sets of 2^eta)
		return [int(math.ceil( n * self.domain_size / (2.0**self.localization_number) )) for n in self.other_oracles]

	def coset_hash_queries(self):
		#return the number of queries for the other oracles
		if not self.check_empty_entries():
			raise CompatibilityError("Evaluating queries for incomplete FRI parameters")

		# return a list whith the same length of other_oracles
		# that in each entry has the nnumber of estimated queries to each 
		# outer oracle

		N = self.domain_size / 2**self.localization_number
		estimated_queries = self.estimate_queries(N)
		return [estimated_queries] * len(self.other_oracles)

	def complete(self):
		if self.query_bound == None or self.domain_dim == None or self.localization_number == None:
			raise CompatibilityError("Unable to deduce FRI parameters")
		else:
			# update the domain size
			self.domain_size = 2**self.domain_dim

			# set the rate of the RS code (hard-coded Aurora's parameters...)
			self.rate = (2*max(self.variables, self.constraints) + 2*self.query_bound)/(2.0**self.domain_dim)
			
			if self.snd_type == "proven":
				# relative distance tested, hard-coded as a function of the rate
				self.relative_distance = min( (1 - 2*self.rate)/2.0, (1 - self.rate)/3.0 )

			elif self.snd_type == "heuristic":
				# See libiop, ldt_reducer, lines 39-42. In this case since the test is ZK
				# max tested rate equals the rate
				self.relative_distance = 1 - self.rate

			# other term in the FRI soundness
			#		d_0 = (1 - 3rho - 2^eta/v)/4 
			# where rho is the rate, eta the localization parameter and v is
			# the domain size
			d_0 = (1.0 - 3.0*self.rate - (2.0**self.localization_number)/(self.domain_size**(0.5)))/4.0

			if self.snd_type == "proven":
				# set the query repetition as the lowest number that achieve the right soundness.
				# soundness is derived from BBHR18. If for the current parameter the soundness bound is
				# 1 or greater we set query_repetition to infinity
				if d_0 <= 0:
					self.query_repetition = float('inf')
				else:
					self.query_repetition = -(self.query_soundness_error)/(np.log2(1 - min(self.relative_distance, d_0)))			
					self.query_repetition = int(math.ceil(self.query_repetition))
			
			elif self.snd_type == "heuristic":
				# set the query repetition as the lowest number that achieve the right heuristic soundness
				# that is (1 - d)^l where d is the relative distance and l is the number of query repetitions.
				# remark that this is always smaller than 1
				self.query_repetition = -self.query_soundness_error / np.log2( 1 - self.relative_distance )
				self.query_repetition = int(math.ceil(self.query_repetition))

	def base_values_from_query_bound(self):
		# Helper for optimize. Set all the independent parameters (beside the query bound) to the minimum value

		if self.query_bound == None:
			raise CompatibilityError("Unable to Reset FRI parameters, query bound set to None")
		else:
			# v > 4 max(n,m) + 4 b because rho > 1/2, otherwise the minimum distance is a negative value
			# with v the domain size, n the number of variables, m the number of constraints, b the queries
			# bound, rho the rate.
			domain_dim = np.log2(4*max(self.variables, self.constraints) + 4*self.query_bound + 1)
			domain_dim = int(math.ceil(domain_dim))

			localization_number = 1
			return(domain_dim, localization_number)

	def estimate_cost(self, estimate = True, very_verbose = False):
		if not self.check_empty_entries():
			raise CompatibilityError("Cannot evalute cost as some fields are set to None")
		else:
			#maximum numer of rounds (BBHR18), recall that rate < 1
			rounds = (self.domain_dim + np.log2(self.rate))/ self.localization_number
			rounds = int(math.ceil(rounds))

			# We initialise the alphabet size, the oracle lenght and query complexity
			# adding at first the masking term. I has
			#	alphabet_size 	= |F| 2^\eta
			#	oracle_length 	= |L| / 2^\eta
			#	queries 		= l
			# with |L| = domain_size, eta = localization_number, l = repetition number
			# |F| the field size

			N = self.domain_size / 2**(self.localization_number)
			estimated_queries = self.estimate_queries(N)

			vec_alphabet_size 	= [self.field_dim * 2**self.localization_number]
			vec_oracle_length 	= [N] 
			vec_query 			= [estimated_queries]			

			for i in range(rounds):
				#number of points in L^(i+1)
				N = self.domain_size / 2**(self.localization_number * (i + 1))
				#number of queries made to f_		
				estimated_queries = self.estimate_queries(N)

				#OPTIMISATION:
				#	if the ratio Q/N is higher than rho then it's more efficient to just
				#	send the polynomial. In that case we truncate FRI recursion
				if estimated_queries/float(N) > self.rate:
					rounds = i
					break

				# add the values to the lists:
				# 	alphabet : |F| * 2^\eta (coset hashing)
				vec_alphabet_size.append(self.field_dim * (2**self.localization_number))
				#	oracle length : |L^(i + i)| (coset hashin)
				vec_oracle_length.append(N)
				#	query complexity: add the estimated number of queries
				vec_query.append(estimated_queries)

			# We add the other oracles coming from outer protocols, considering coset hashing
			vec_alphabet_size = self.coset_hash_alphabet() + vec_alphabet_size
			vec_oracle_length = self.coset_hash_length() + vec_oracle_length
			vec_query = self.coset_hash_queries() + vec_query

			# For the round complexity, from left to right we have
			# 	the round executed before running the LDT (number of entries in other_oracles)
			# 	one round for the LDT masking term
			# 	the FRI reduction rounds
			# 	the final direct LDT round (after reduction)
			total_rounds = len(self.other_oracles) + 1 + rounds + 1

			bsc_cost = com_BCS(self.hash_size, vec_alphabet_size, 
				total_rounds, vec_oracle_length, vec_query, estimate)
			
			# Dimension of the domain over which we perform the direct LDT
			last_domain_size = self.domain_size / 2**(self.localization_number * rounds)
			# Cost of the last, direct LDT, done by sending all the coefficient of the polynomial
			direct_ldt_costs = self.field_dim * self.rate * last_domain_size
			

			#if very_verbose:
				# notice that with VERBOSE_FLAG = False, nothing is printed anyway
			#	print_cost("BSC in FRI", bsc_cost)
			#	print_cost("Direct LDT in FRI", direct_ldt_costs)
			#print_cost("Total FRI", bsc_cost + direct_ldt_costs)

			return bsc_cost + direct_ldt_costs

	def avg_cost(self):
		# Run estimate_cost setting the estimate flag to False.
		# this cause BSC to run the (expensive) Montecarlo simulation
		return self.estimate_cost(estimate = False, very_verbose = True)

	def optimize(self):
		# idea taken from libiop: we set the query bound to 0 and keep runnning
		# the optimizer until the query bound is strictly smaller than the 
		# query repetitions
	
		self.query_bound = 0
		tuple_out = None # Best triplet of parameters (q, eta, bv)
		cost_out = float('inf') # current best cost
		cost_tmp = float('inf') # current cost
		
		#to avoid repetitions we exclude already tested values
		tested_list = []

		while(True): #Emulate a do-while to find the best condition
			#set the minimum possible values for eta and bv
			domain_dim_0, localization_number_0 = self.base_values_from_query_bound() 
			exit_flag = True

			# we fix an arbitrary bound on the domain dimension to keep the computation
			# under controll.
			for self.domain_dim in range(domain_dim_0, domain_dim_0 + FRI_MAX_DOMAIN_DIM + 1):
				# we make the loc. number run up to a fixed bound 
				for self.localization_number in range(localization_number_0, FRI_MAX_LOCALIZATION_NUM + 1):

					# we compute the remaining parameters
					self.complete()

					if self.query_repetition != float('inf') and \
					   (self.domain_dim, self.localization_number) not in tested_list:

						# if the query bound is not smaller than the number of queries made to the base oracle:
						if self.query_bound >= self.query_repetition * 2**self.localization_number:
							
							tested_list.append((self.domain_dim, self.localization_number))

							cost_tmp = self.estimate_cost()
							if cost_tmp < cost_out:
								cost_out = cost_tmp
								tuple_out = (self.query_bound, self.domain_dim, self.localization_number)
							
							print_verbose_message("\tTrying (b, dim L, eta) = ({:d}, {:d}, {:1d})".format(
									self.query_bound, self.domain_dim, self.localization_number))
							print_verbose_cost("\tFinal proof size", cost_tmp)
							print_verbose_message("")

						else:
							exit_flag = False 
							# there is still some couple (dimL, eta) not captured by the
							#  current estimated queries bound



			# Exit condition:
			# When he manages to match the query bound from the given query_repetitions
			if exit_flag:
				break
			else:
				self.query_bound += 1

		self.query_bound, self.domain_dim, self.localization_number = tuple_out
		self.complete()

	def optimal_cost(self):
		self.optimize()
		return self.avg_cost()

# - - - - - - - - - - - - - - - - - #
#	Aurora's parameters and test 	#
# - - - - - - - - - - - - - - - - - #

class AURORA_parameters:
	def __init__(self, n, m, fd, sp, snd_type, protocol_type, rmfe = None):
		# n : number of input gates
		# m : number of AND gates
		# fd : field dimension (bitsize)
		# sp : security parameter

		# Type of protocol measured:
		#	"standard" : classic Aurora
		#	"optimised" : our version

		self.field_dim = fd
		self.security_parameter = sp

		self.protocol_type = protocol_type 	
		self.snd_type = snd_type

		if self.protocol_type == "standard":
			#set the standard parameters for standard AURORA
			self.RMFE = None
			self.variables = n 			# input gates stays constant
			self.constraints = m*2		# the number of constraints is doubled if using x^2 - x
		
		elif self.protocol_type == "optimised":
			#set the parameters for our optimised version
			self.RMFE = rmfe[0]												# rmfe = (RMFE, field_dim)
			self.field_dim = rmfe[1]										#  ...
			self.variables = 2**math.ceil(np.log2( float(n)/self.RMFE ))	# n' = n/k approximated to
			self.variables = int(self.variables)							#  the smallest pow of 2
			self.constraints = 2**math.ceil(np.log2( float(m)/self.RMFE ))	# m' = m/k approximated to
			self.constraints = int(self.constraints)						#  the smallest pow of 2

		self.FRI_parameters = FRI_parameters(\
			self.variables,					# Variable number (n)					\
			self.constraints,				# Constraints number (m)				\
			self.field_dim,					# Field dimension (fd)					\
			self.security_parameter + 3,	# Interactive Soundness Error (isp)		\
			self.security_parameter + 1, 	# Query Soundness Error (qsp)			\
			self.security_parameter*2, 		# Hash size (h)							\
			snd_type = snd_type, other_oracles = self.oracles_number())

	def __str__(self):
		return str(self.FRI_parameters)

	def oracles_number(self):
		# Number of oracles sent during the protocol. In this we consider
		# 	round 1: f_w, f_Az, f_Bz, f_Cz, sumcheck_mask
		#	round 2: sumcheck_rest
		if self.protocol_type == "standard":
			return [5, 1]
		elif self.protocol_type == "optimised":
			return [6, 1]

	def extra_communication(self):
		# Ammount of data sent by the prover before the LDT and NOT
		#  in the form of an oracle. In theory this consist of one 
		#  element in the sumcheck protocol, but this element can
		#  be avoided is the prover already knows the sum AND the
		#  masking term is chosen to have zero sum. In the optimised version
		#  this consist of the vectors in the Modular Lincheck
		if self.protocol_type == "standard":
			return 0
		elif self.protocol_type == "optimised":
			return 2 * self.field_dim * (self.security_parameter + 3)

	def optimal_cost(self):
		# Optimal cost for Aurora with given parameters
		fri_cost = self.FRI_parameters.optimal_cost()
		extra_cost = self.extra_communication()
		return fri_cost + extra_cost

# - - - - - - - - - - - - - - - #
# 	Testing functionalities 	#
# - - - - - - - - - - - - - - - #

def _parameters_class(scheme):
	# Helper function for comparison_to_csv
	#  return the parameter class requested from the test
	#  
	#  scheme : "ligero" or "aurora"

	if scheme == "ligero":
		return LIGERO_parameters
	elif scheme == "aurora":
		return AURORA_parameters
	else:
		# debug
		assert False, "Requested unkown parameters class"

def comparison_to_csv(scheme, filename, dim_min, dim_max, fd, rmfe, sp, snd_type):
	# Make a comparison between plain scheme over F2 and our optimised version
	#
	#  scheme 		: "aurora" or "ligero" - choose which scheme to test
	#  filename 	: name of the file in which final data is stored
	#  dim_min 		: input size in the test start from 2**dim_min
	#  dim_max 		: input size in the test end with 2**dim_max (inclusive)
	#  fd 			: Aurora's field dimension
	#  rmfe 		: A tuple (k,m) of rmfe parameters used in the opt
	#  sp 			: security parameter
	#  snd_type 	: soundness type

	with open(filename + ".csv", "w") as file:
		# list of the constraints tested for both the standard and optimised version of scheme
		vec_constraints_dim = [n for n in range(dim_min, dim_max + 1)]
		
		print_message("Writing results in {}\n".format(filename + ".csv"))
		vec_to_csv(file, vec_constraints_dim)

		# choose the parameter class
		parameters_class = _parameters_class(scheme)

		vec_results = []
		for n in vec_constraints_dim:
			# execute the standard test
			print_message("Running standard {:s} for n = {:2d}".format(scheme, n))
			parameters = parameters_class(2**n, 2**n, fd, sp, snd_type, "standard")
			vec_results.append(parameters.optimal_cost())

			# if VERY_VERBOSE_FLAG
			print_verbose_cost("\nStandard {:s} with {:s} soundness, n = {:2d}".format(scheme, snd_type, n), \
				vec_results[-1])
			print_separator(very_verbose = True)

		# store the results
		print_verbose_message("End of standard {:s} tests".format(scheme))
		print_separator()
		vec_to_csv(file, vec_results)

		vec_results = []
		for n in vec_constraints_dim:
			# execute the optimised test
			print_message("Running optimised {:s} for n = {:2d}, rmfe = {}".format(scheme, n, rmfe))
			parameters = parameters_class(2**n, 2**n, fd, sp, snd_type, "optimised", rmfe = rmfe)
			vec_results.append(parameters.optimal_cost())

			# if VERY_VERBOSE_FLAG
			print_verbose_cost("\nOptimised {:s} with {:s} soundness, n = {:2d}".format(scheme, snd_type, n), \
				vec_results[-1])
			print_separator(very_verbose = True)

		# store the results
		print_verbose_message("End of optimised {:s} tests".format(scheme))
		vec_to_csv(file, vec_results)

def AURORA_test(dim = 18, fd = 192, sp = 128, rmfe = None, snd_type = "proven", protocol_type = "standard"):
	print_message("testing {:s} aurora, n = {:2d}".format(protocol_type, dim))
	aurora_p = AURORA_parameters(2**dim, 2**dim, fd, sp, rmfe = rmfe, snd_type = snd_type, protocol_type = protocol_type)

	cost = aurora_p.optimal_cost()
	print_message(str(aurora_p))
	return cost

#AURORA_comparison_to_csv("test_aurora_heuristic_rmfe_48_160", snd_type = "heuristic", rmfe = (48, 160))

#LIGERO_comparison_to_csv("Data/ligero_proven_RMFE_45_145", rmfe = [(45, 145)])

if __name__ == "__main__":
	argv = sys.argv[1:]
	
	# -l 	--ligero 				Runs a comparison for Ligero
	# -a 	--aurora 				Runs a comparison for Aurora
	# -tl   --test_ligero			Runs a test for Ligero
	# -ta 	--test_aurora 			Runs a test for Aurora
	#
	# -sp 	--security_parameter	Set the security parameter (1 arg)
	# -ld 	--lowest_dimension		Set the minimum dimension tested, inclusive (1 arg)
	# -hd 	--highest_dimension 	Set the maximum dimension tested, inclusive (1 arg)
	# -h 	--heuristic 			Use heuristic soundness bounds (0 arg)
	# -p 	--proven 				Use proven soudness bounds (0 arg)
	# -r 	--rmfe 					Set the RMFE for the optmised version (2 arg)
	# -fd 	--field_dimension 		Set the field dimension for the standard version (1 arg)
	# -f 	--file_name 			Choose the file name (1 arg)
	# -D 	--directory				Absolute or relative path, added to filename as a prefix (1 arg)
	# -v 	--verbose 				Print some information on the screen
	# -vv 	--very_verbose 			Print more informations on the screen
	#
	# maybe make a parameter class to handle this in a cleaner way
	# clearly all but verbose flags

	test_type = "aurora"
	sp = 128
	dim_min = 8
	dim_max = 20
	rmfe = None 		#auto completed after parsing
	snd_type = "proven" 
	fd = None			#auto completed after parsing
	filename = None 	#auto completed after parsing
	dirname = ""				
	VERBOSE_FLAG = False
	VERY_VERBOSE_FLAG = False
	DEBUG = False 

	# look for the help flag
	if "-help" in argv:
		print(str_help())
		sys.exit(0)

	# Parse the input and set the variables
	try:
		n, i = len(argv), 0
		# check the arguments received
		while i < n:
			key = argv[i]
			
			if key in ["-l", "--ligero"]:
				test_type = "ligero"
			
			elif key in ["-a", "--aurora"]:
				test_type = "aurora"

			elif key in ["-tl", "--test_ligero"]:
				test_type = "test_ligero"
				VERBOSE_FLAG = True

			elif key in ["-ta", "--test_aurora"]:
				test_type = "test_aurora"
				VERBOSE_FLAG = True

			elif key in ["-sp", "--security_parameter"]:
				sp = int(argv[i+1])
				i += 1

			elif key in ["-ld", "--lowest_dimension"]:
				dim_min = int(argv[i+1])
				i += 1

			elif key in ["-hd", "--highest_dimension"]:
				dim_max = int(argv[i+1])
				i += 1

			elif key in ["-h", "--heuristic"]:
				snd_type = "heuristic"

			elif key in ["-p", "--proven"]:
				snd_type = "proven"

			elif key in ["-r", "--rmfe"]:
				rmfe = ( int(argv[i+1]), int(argv[i+2]) )
				i += 2

			elif key in ["-fd", "--field_dimension"]:
				fd = int(argv[i+1])
				i += 1

			elif key in ["-f", "--file_name"]:
				filename = str(argv[i+1])
				i += 1

			elif key in ["-D", "--directory"]:
				dirname = str(argv[i+1])
				i += 1

			elif key in ["-v", "--verbose"]:
				VERBOSE_FLAG = True

			elif key in ["-vv", "--very_verbose"]:
				VERY_VERBOSE_FLAG = True
				VERBOSE_FLAG = True

			else:
				# If there is no flag where a flag was expected
				#  raise an InputError 
				raise InputError

			i += 1
	except (IndexError, ValueError, InputError):
		print(str_synopsis())
		sys.exit(0)

	# set rmfe
	if rmfe == None:
		if test_type in ["aurora", "test_aurora"]:
			rmfe = (48, 198)
		elif test_type in ["ligero", "test_ligero"]:
			rmfe = (48, 160)

	# set fd
	if test_type in ["aurora", "test_aurora"]:
		fd = 192

	# set filename
	if filename == None:
		filename = "{:s}_{:s}_sp_{:d}_rmfe_{:d}-{:d}".format(test_type, snd_type, sp, *rmfe)

	# test execution
	if test_type in ["aurora", "ligero"]:
		comparison_to_csv(test_type, filename, dim_min = dim_min, dim_max = dim_max, fd = fd, rmfe = rmfe, \
			sp = sp, snd_type = snd_type)

	elif test_type == "test_aurora":
		for d in range(dim_min, dim_max + 1):
			cost = AURORA_test(dim = d, fd = fd, sp = sp, snd_type = snd_type, protocol_type = "standard")
			print_cost("Final cost", cost)
			print_separator()
		for d in range(dim_min, dim_max + 1):
			cost = AURORA_test(dim = d, fd = fd, sp = sp, rmfe = rmfe, snd_type = snd_type, protocol_type = "optimised")
			print_cost("Final cost", cost)
			print_separator()
		print_message("End of test")

	elif test_type == "test_ligero":
		for d in range(dim_min, dim_max + 1):
			cost = LIGERO_test(dim = d, fd = fd, sp = sp, snd_type = snd_type, protocol_type = "standard")
			print_cost("Final cost", cost)
			print_separator()
		for d in range(dim_min, dim_max + 1):
			cost = LIGERO_test(dim = d, fd = fd, sp = sp, rmfe = rmfe, snd_type = snd_type, protocol_type = "optimised")
			print_cost("Final cost", cost)
			print_separator()
		print_message("End of test")