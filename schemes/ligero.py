from .constants import *
from .bcs import *
from ._decorators import *

from util.utilities import *

# - - - - - - - - - - - - - - - #
#	Ligero for binary circuits	#
# - - - - - - - - - - - - - - - #

@protocol_types("booligero")
class ligero_parameters:
	def __init__(self, n, m, fd, sp, snd_type, protocol_type, rmfe = None):

		self.protocol_type = protocol_type	# "standard"/"optimised"/"booligero"
		self.field_dim = fd 				# log2(|F|)	

		if self.protocol_type == "optimised":
			# debug - the rmfe parameter has to be given
			assert rmfe != None, "LIGERO, rmfe not given in optimised version"
			# - - - - - end of debug - - - - - #
			self.rmfe, self.field_dim = rmfe

		# checking if the size of the field is fixed. A not given
		#  value implies that the field size can be adjusted to minimize
		#  communication costs.
		if self.field_dim == None:
			self.field_type = "variable"
			self.field_size = None						
		else:
			self.field_type = "fixed"
			self.field_size = 2**self.field_dim

		# Protocol-dependent parameters
		if self.protocol_type == "standard":
			
			# For each variable we need the constraint x^2 = x
			self.variables = n
			self.constraints = m + n
			
			# Overall soundness is the sum of isp and qsp 
			self.interactive_soundness_error = sp + 1
			self.query_soundness_error = sp + 1

			# Range in which we look for optimal l
			self.LGR_MIN_L = 1/10.0
			self.LGR_MAX_L = 2/3.0

		elif self.protocol_type == "optimised":

			# Variables and constrains scale down of a rmfe factor
			self.variables = ceil(n/self.rmfe)
			self.constraints = ceil(m/self.rmfe)

			# Overall soundness is the sum of 3 isp and one qsp
			self.interactive_soundness_error = sp + 3
			self.query_soundness_error = sp + 1

			# Range in which we look for optimal l
			self.LGR_MIN_L = 1/6.0
			self.LGR_MAX_L = 2.5

		elif self.protocol_type == "booligero":

			# Variables and constraints are recorded here just for __str__
			self.variables = n
			self.constraints = m

			def w0(f):
				return int(math.floor(f**(0.5)))

			def w1(f = None):
				if f is None:
					f = self.field_dim
				return ceil(f/self.w0(f))

			self.w0 = w0
			self.w1 = w1

			# Overall soundness is the sum of 2 ips and one qsp
			self.interactive_soundness_error = sp + 2
			self.query_soundness_error = sp + 1

			# Range in which we look for optimal l
			self.LGR_MIN_L = 1/12.0
			self.LGR_MAX_L = 1.8

		# Other parameter fixed ahead
		self.hash_size = 2*sp
		self.security_parameter = sp 	 					# kappa
		self.snd_type = snd_type

		# Parameters that look for in order to find optimal costs
		self.domain_dim = None				# n
		self.domain_size = None 			# 2^n 	**
		self.degree = None					# k 		
		self.proximity_parameter = None 	# e 	**
		self.queries = None 				# t 	**
		self.factor_l = None				# l 		
		self.factor_m1 = None				# m1 	** codewords to fill w
		self.factor_m2 = None				# m2 	** cws to fill x, y, z, [t]
		self.interactive_repetitions = None	# sigma **
			

	def check_empty_entries(self, var = None, exclude = False):
		# check that all the variables in var are not set to None
		#  if exclude then check all variables beside those
		#  whose name is in var
		# 
		# It returns the name of the empty variable. If no empty
		#  variable is found, it returns None
		if var is None:
			var = []

		for key, value in vars(self).items():
			if value == None and xor(key in var, exclude):
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
			out += "RMFE parameters:\t({}, {})\n".format(self.rmfe, self.field_dim)

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

	def compute_queries(self):
		'''Returns the number of queries required to have a query soundness
		error smaller that 2^{-qsp}
		'''

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
			sigma = lambda f : ceil( isp/(f - pp) )
		
		elif self.snd_type == "heuristic":
			# assuming an interactive soundness error of (e + 1)/|F|^sigma
			sigma = lambda f : ceil( (isp + pp)/f )	
		
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
		for f in range(self.domain_dim, ceil(isp + pp + 1)):

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

		if self.protocol_type == "standard":
			# We should set 3*(k + l - 1) as we only perform 3 lincheck in the
			#  R1CS, however ligero allows for linear test to be batched
			#  without any soundness loss
			#
			# Costs ammount to
			#
			# LDT :	n 			: a codeword)
			# Lin : k + l - 1 	: a k + l - 2 degree polynomial
			# Row : 2k - 1 		: a 2k - 2 degree polynomial
			out = f*sigma*(n + (k + l - 1) + (2*k - 1))

		if self.protocol_type == "optimised":
			# On top of the standard cost for performing LDT, Lin and Row we
			#  also send 2 vector whose length is roughly the size of the s.p.
			out = f*sigma*(n + (k + l - 1) + (2*k - 1))
			out += 2*f*self.interactive_soundness_error

		if self.protocol_type == "booligero":
			# On top of the standard cost for performing LDT, Lin and Row we
			#  also send 3 vectors to check the patterns of size
			#  f * (w1 + 1) * sp
			out = f*sigma*(n + (k + l - 1) + (2*k - 1))
			out += 3*f*(self.w1(f) + 1)*self.interactive_soundness_error

		if query_flag:
			out += self.queries * self.find_alphabet_size(
				sigma = sigma, f = f, sum_flag = True)

		return out

	def find_rounds(self):
		if self.protocol_type == "standard":
			return 2

		elif self.protocol_type == "optimised":
			return 3

		elif self.protocol_type == "booligero":
			return 3

	def find_alphabet_size(self, sigma = None, f = None, sum_flag = False):
		''' return the alphabet size of the (column compressed) oracles 
		used in BSC.
		
		sigma and f are provided to make "speculative" calls. if not 
		provided we set them as the associated values in self

		if sum_flag is True returns the sum of the sizes of all the oracles sent
		(if the prover sends more than one oracle)
		'''

		# TODO - remove in the next commit
		# debug - check for empty entries
		assert self.check_empty_entries(var = ("protocol_type", "factor_m1", "factor_m2")) == None, \
			"LIGERO, empty required variable"

		if sigma == None:
			# TODO - remove in the next commit
			# debug - check for empty entries
			assert self.check_empty_entries(var = ("interactive_repetitions")) == None, \
				"LIGERO, empty required variable"
			# - - - - - end of debug - - - - - #
			sigma = self.interactive_repetitions

		if f == None:
			# TODO - remove in the next commit
			# debug - check for empty entries
			assert self.check_empty_entries(var = ("field_dim")) == None, \
				"LIGERO, empty required variable"
			# - - - - end of debug - - - - - #
			f = self.field_dim

		##
		if self.protocol_type == "standard":
			# Oracles: w (m1), x, y, z (m2)
			# Masking terms: sigma for 1 LDT, 1 Lin, 1 Row
			out = (self.factor_m1 + 3*self.factor_m2 + 3*sigma)*f

		elif self.protocol_type == "optimised":
			# Oracles: w (m1), x, y, z, t (m2)
			# Masking terms: sigma for 1 LDT, 1 Lin_h, 1 Row
			out = (self.factor_m1 + 4*self.factor_m2 + 3*sigma)*f

			# TODO - remove with the next commit
			# debug - check for empty entries
			assert self.check_empty_entries(var = ("interactive_soundness_error", "factor_l")) == None, \
				"LIGERO, empty required variable"
			# - - - - - end of debug - - - - - #

			# masking terms for the modular lincheck
			out += ceil((3*self.interactive_soundness_error)/self.factor_l)*f

		elif self.protocol_type == "booligero":
			# Oracles
			#  Let n = var, m = cont, f = field dimension, sp = security param
			#	
			# Round 1
			#  w 			: n/f*l 			f.e.
			#  x, y, z 		: m/f*l 			f.e.
			#  x^, y^, z^	: (m/f*l)*w1 		f.e.
			#  mask 		: 3 sigma  			f.e.
			# 
			# Round 2 							
			#  ax, ay, az 	: 3*(w1 + 1)*sp 	f.e.

			# Number of field elements on a single caracter
			round1 = ceil((self.variables//f
				+ 3*self.constraints//f
				+ 3*self.w1(f)*self.constraints//f) / self.factor_l)
			# Add the masking terms
			round1 += 3*sigma
			# Convert number of field elements to number of bits
			round1 *= f

			# Number of bits on a signle caracter
			round2 = 3*(self.w1(f) + 1)*f*self.interactive_soundness_error

			out = [round1, round2]

			# Messy...
			if sum_flag:
				return sum(out)


		return out

	def find_oracle_length(self):
		if self.protocol_type == "standard":
			return self.domain_size

		elif self.protocol_type == "booligero":
			return [self.domain_size, self.domain_size]

		elif self.protocol_type == "optimised":
			return self.domain_size

	def find_queries(self):
		if self.protocol_type == "standard":
			return self.queries

		elif self.protocol_type == "booligero":
			return [self.queries, self.queries]

		elif self.protocol_type == "optimised":
			return self.queries

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
		self.factor_m1 = self.variables / self.factor_l
		self.factor_m1 = ceil(self.factor_m1)

		# m2 = ceil(m/l) so that l m2 >= n
		self.factor_m2 = self.constraints / self.factor_l
		self.factor_m2 = ceil(self.factor_m2)

	def complete(self):
		# Using degree, domain_dim and factor_l set the other parameters
		#  according to [AHIV17]

		# debug - check for empty variables
		assert self.check_empty_entries(var = ("domain_dim", "factor_l", "degree")) == None, \
			"LIGERO empty required variable"
		# - - - - - - - - - - end of debug - - - - - - - - - - #

		# number of points in the RS domain
		self.domain_size = 2**self.domain_dim

		# e = floor((n - 2k + 1)/3) with strong Ligero
		#
		self.proximity_parameter = (self.domain_size - 2*self.degree + 1)/3.0
		self.proximity_parameter = int(math.floor(self.proximity_parameter))

		# set the factors m1 and m2
		self.set_factor_m()

		# t is the minimum value that makes the query soundness error smaller
		#  than 2^(-qsp)
		self.queries = self.compute_queries()
		
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
				# Sigma is the minimum value that makes the interactive
				#  soundness error smaller than 2^(-isp)
				self.interactive_repetitions = (self.interactive_soundness_error
					/ (self.field_dim - np.log2(self.proximity_parameter + 1)))
				self.interactive_repetitions = ceil(
					self.interactive_repetitions)
			
			elif self.snd_type == "heuristic":
				self.interactive_repetitions = (self.interactive_soundness_error
					+ np.log2(self.proximity_parameter + 1))/self.field_dim
				self.interactive_repetitions = ceil(
					self.interactive_repetitions)

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
		oracle_length = self.find_oracle_length()
		queries = self.find_queries()

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
		l_min = ceil(l_0 * self.LGR_MIN_L)
		l_max = ceil(l_0 * self.LGR_MAX_L)
		l_step = ceil((l_max - l_min)/LGR_MAX_SAMPLE)

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
			domain_values = [ceil(np.log2(2*self.degree + 1)) + i 
				for i in range(LGR_MAX_DOMAIN_DIM)]
			domain_values_tmp = []

			while(domain_values):								# while there are still values to check
				min_degree_bound = float("inf")

				for self.domain_dim in domain_values:

					self.complete()								# set the other parameters as functions of m,k,n
					if self.degree - self.factor_l >= self.queries:
						print_verbose_message("\tTesting (l, k, log2n, delta) = ({}, {}, {:2d}, {})".format(
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
