from .constants import *
from .fri import *

from util import *

# - - - - - - - - - - - - - - - - - #
#	Ligero++'s parameters and test	#
# - - - - - - - - - - - - - - - - - #

class ligeropp_parameters:
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

		self.query_soundness_error = sp
		self.interactive_soundness_error = sp

		if self.protocol_type == "standard":
			# Set the standard parameters for standard ligero++
			self.RMFE = None
			# var = n, constraints = m + n (for boolean circuit we need to add
			#  x^2 - x for each variable
			self.variables = n 			
			self.constraints = m + n
		
		elif self.protocol_type == "optimised":
			# Set the parameters for our optimised version
			#  rmfe = (k, fd)
			self.RMFE, self.field_dim = rmfe

			# var = n/k, con = m/k
			self.variables = ceil(n/self.rmfe)
			self.constraints = ceil(m/self.rmfe)

		self.fri_parameters = fri_parameters(
			lambda q : max(self.factor_m1, 2*self.factor_m2), # deg_sigma
			lambda q : 2*max(self.factor_m1, self.factor_m2), # deg_rho
			self.field_dim,					# Field dimension (fd)
			self.security_parameter + 3,	# Interactive Soundness Error (isp)
			self.security_parameter + 1, 	# Query Soundness Error (qsp)
			self.security_parameter*2, 		# Hash size (h)
			snd_type=snd_type,
			zk=False)						# Zero-Knowledge flag [not required]

		# Free parameters
		self.domain_dim = None				# log n ###
		self.domain_size = None				# n
		self.degree = None					# k
		self.factor_log_l = None			# log l ###
		self.factor_l = None				# l
		self.factor_m1 = None				# m1
		self.factor_m2 = None				# m2
		#self.RS_min_distance = None		# d
		#self.proximity_parameter = None	# e
		self.queries = None					# t

		# In the paper two version are presented
		#  "zk-ligero" uses techniques from ligero to achieve zero knowledge
		#  "zk-IPA" uses a zk low degree test and a zk sumcheck to win
		self.mode = None

	def __str__(self):
		out = "\n- - - - L++ pramaters - - - -\n"
		out += "Variables:\t\t{}\n".format(self.variables)
		out += "Constraints:\t\t{}\n".format(self.constraints)
		out += "Domain dimension:\t{}\n".format(self.domain_dim)
		out += "Degree used:\t\t{}\n".format(self.degree)
		out += "Factor l:\t\t{}\n".format(self.factor_l)
		out += "Factor m1:\t\t{}\n".format(self.factor_m1)
		out += "Factor m2:\t\t{}\n".format(self.factor_m2)
		out += "IPAs queried:\t\t{}\n".format(self.queries)

		out += str(self.fri_parameters)

		return out

	def oracles_number(self):
		# Number of oracles sent during the protocol. In this case we consider
		# 	round 1: f_w (l), f_x (l), f_Bz (l), f_Cz (l)
		#
		# Moreover in mode = zk-ligero
		# 	round 1: sumcheck_rest
		#
		# In mode = zk-IPA
		#	round 1: sumcheck_mask
		#	round 2: sumcheck_rest
		#
		#  Remember that the sumcheck does not need to be zk - hence only the
		#  remainder is sent by the prover

		if self.protocol_type == "standard":
			if self.mode == "zk-ligero":
				return [4*self.queries + 1]
			elif self.mode == "zk-IPA":
				return [4*self.queries + 1, 1]

		elif self.protocol_type == "optimised":
			assert False, "Value not yet defined"
			#return [5*self.queries + 1]

	def extra_communication(self):
		# Ammount of data sent by the prover before the LDT and NOT
		#  in the form of an oracle. In theory this consist of one 
		#  element in the sumcheck protocol, but this element can
		#  be avoided is the prover already knows the sum AND the
		#  masking term is chosen to have zero sum. In the optimised version
		#  this consist of the vectors in the Modular Lincheck

		if self.protocol_type == "standard":
			n = self.domain_size
			k = self.degree
			l = self.factor_l
			fd = self.field_dim
			return fd*(n + 3*(k + l - 1) + (2*k - 1))

		elif self.protocol_type == "optimised":
			assert False, "Value not yet defined"
			#return 2 * self.field_dim * (self.security_parameter + 3)

	def _query_soundness_error(self, n, l, t):
		'''Return the soundness error given
		n : domain size
		l : factor_l
		t : number of queries
		'''

		#DEBUG - old version
		#err_ldt = ((2*n + l + t - 1)/(3*n))**t
		#err_row = ((n + 5*l + 5*t + 1)/(3*n))**t

		k = l + t
		e = (n - k + 1)//3

		err = max(combr(n - e - 1, n, t), combr(e + 2*k - 2, n, t))
		return err

	def complete(self):
		'''if self.domain_dim and self.factor_log_l are defined, compute
		the remaining attributes
		'''

		self.domain_size = 2**self.domain_dim
		self.factor_l = 2**self.factor_log_l

		# we evaluate by binary search [from util.utilities]
		err = 2**(-self.query_soundness_error)
		n = self.domain_size
		l = self.factor_l
		fun = lambda t : self._query_soundness_error(n, l, t) <= err
		
		self.queries = binsearch(fun, 1, upper_bound = n - 1)

		if self.queries == float('inf'):

			# if we would need more queries than possibile in zk-ligero, abort
			if self.mode == "zk-ligero":
				 return

			# if in zk-IPA, query every column
			elif self.mode == "zk-IPA":
				self.queries = self.domain_size

		# Next compute k
		if self.mode == "zk-ligero":
			self.degree = self.factor_l + self.queries
		else:
			self.degree = self.factor_l

		# compute m1, m2
		self.factor_m1 = ceil(self.variables/self.factor_l)
		self.factor_m2 = ceil(self.constraints/self.factor_l)

		# Finally, update FRI's parameters
		self.fri_parameters.other_oracles = self.oracles_number()
		self.fri_parameters.optimize()
		if self.mode == "zk-ligero":
			self.fri_parameters.zero_knowledge = False
		elif self.mode == "zk-IPA":
			self.fri_parameters.zero_knowledge = True

	def estimate_cost(self, estimate = True):
		'''Returns the cost of applying Ligero++ with current parameters
	
		If estimate is True the cost is computed faster but in a approximate way
		Otherwise the computation is slower but the result more accurate
		'''
		fri_cost = self.fri_parameters.estimate_cost(estimate = estimate)
		extra_cost = self.extra_communication()

		#DEBUG
		#print_cost("fri's cost", fri_cost)
		#print_cost("extra cost", extra_cost)
		#print("")
		#DEBUG

		return fri_cost + extra_cost

	def avg_cost(self):
		return self.estimate_cost(estimate = False)

	def optimize(self):

		# We perform ax exhaustive search for domain_dim and log2_l
		#  Observe that n >= 2k = 2l + 2t. Hence we have that
		#
		#  .:	log(n) > log(2l) = log(l) + 1
		#  ->	log(n) >= log(l) + 2
		#

		tuple_out = None
		cost_tmp = float('inf')
		cost_out = float('inf')

		tuples = [(mode, log_n, log_l) 
			for mode in ("zk-ligero", "zk-IPA")
			for log_l in range(1, LPP_MAX_LOG_L)
			for log_n in range(log_l + 2, log_l + 2 + LPP_MAX_DOMAIN_DIM)
		]

		for (mode, log_n, log_l) in tuples:
			self.mode = mode
			self.domain_dim = log_n
			self.factor_log_l = log_l
			self.complete()

			# We filter out tuples that cannot achieve the right soundness error
			if self.queries == float('inf'):
				continue

			cost_tmp = self.estimate_cost(estimate = True)

			if cost_out > cost_tmp:
				cost_out = cost_tmp
				tuple_out = (mode, log_n, log_l)

				#DEBUG
			print_message("\tTested tuple: {}".format((mode, log_n, log_l)))
			print_cost("\tCost", cost_tmp)
			print_message("")

			print_verbose_message(
				"\tTrying (dim L, dim H) = ({:2d}, {:2d})".format(log_n, log_l))
			print_verbose_cost("\tFinal proof size", cost_tmp)
			print_verbose_message("")

		# set the optimal value found so far
		self.mode, self.domain_dim, self.factor_log_l = tuple_out
		self.complete()

	def optimal_cost(self):
		self.optimize()
		# DEBUG
		return self.estimate_cost(estimate = False)








