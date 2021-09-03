from .constants import *
from .bcs import *

from util.utilities import *

# - - - - - - - - - - - - - - - - - #
#	FRI soundness error and costs	#
# - - - - - - - - - - - - - - - - - #

class fri_parameters:
	def __init__(self, n, m, fd, isp, qsp, h, snd_type = "proven", 
		other_oracles = None):

		if other_oracles is None:
			other_oracles = []

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

	def check_empty_entries(self, var = None, exclude = False):
		# check if any field is empty - should be a method in a super class
		#
		if var is None:
			var = []

		for key, value in vars(self).items():
			if value == None and xor(key in var, exclude):
				return key
		else:
			return None

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

		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["query_repetition"]) == None, "FRI, missing used variable"
		# - - - end of debug - - - #

		estimated_queries = 0.0

		for j in range(FRI_ITER):									
			sample = rnd.randint(N, size = self.query_repetition)	#sample l points
			estimated_queries += len(set(sample))					#consider only the distinct ones

		estimated_queries = int(math.ceil(estimated_queries/float(FRI_ITER)))
		return estimated_queries

	def coset_hash_alphabet(self):
		# Set the alphabet size for the other oracles (leaves are grouped in sets of 2^eta)

		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["field_dim", "localization_number", "other_oracles"]) == None, \
			"FRI, missing required variable"
		# - - - - end of debug - - - - #
		return [int(math.ceil( n * self.field_dim * (2.0**self.localization_number) )) for n in self.other_oracles]

	def coset_hash_length(self):
		# Set the proof length for the other oracles (leaves are grouped in sets of 2^eta)

		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["domain_size", "localization_number", "other_oracles"]) == None, \
			"FRI, missing required variable"
		# - - - - end of debug - - - - #

		return [int(math.ceil( n * self.domain_size / (2.0**self.localization_number) )) for n in self.other_oracles]

	def coset_hash_queries(self):
		#return the number of queries for the other oracles
		
		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["domain_size", "localization_number", "other_oracles"]) == None, \
			"FRI, missing required variable"
		# - - - - end of debug - - - - #

		# return a list whith the same length of other_oracles
		# that in each entry has the nnumber of estimated queries to each 
		# outer oracle

		N = self.domain_size / 2**self.localization_number
		estimated_queries = self.estimate_queries(N)
		return [estimated_queries] * len(self.other_oracles)

	def complete(self):
		# TODO ... add a comment ...

		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["query_bound", "domain_dim", "localization_number"]) == None, \
			"FRI, missing required variable"
		# - - - - end of debug - - - - #
			
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

		# - - - - - debug - - - - - #
		assert self.check_empty_entries(var = ["variables", "constraints", "query_bound"]) == None, \
			"FRI, missing required variable"
		# - - - - end of debug - - - - #		

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
		# - - - - - debug - - - - - #
		assert self.check_empty_entries(exclude = True) == None, "FRI, missing required variables"
		# - - - - end of debug - - - - #

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

		return bsc_cost + direct_ldt_costs

	def avg_cost(self):
		# Run estimate_cost setting the estimate flag to False.
		# this cause BSC to run the (expensive) Montecarlo simulation
		return self.estimate_cost(estimate = False, very_verbose = True)

	def optimize(self):
		# idea taken from libiop: we set the query bound to 0 and keep runnning
		#  the optimizer until the query bound is smaller than the 
		#  query repetitions 
	
		self.query_bound = 0
		tuple_out = None 			# Best triplet of parameters (q, eta, bv)
		cost_out = float('inf') 	# current best cost
		cost_tmp = float('inf') 	# current cost
		
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
