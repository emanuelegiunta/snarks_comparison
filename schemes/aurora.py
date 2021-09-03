from .constants import *
from .fri import *

# - - - - - - - - - - - - - - - - - #
#	Aurora's parameters and test 	#
# - - - - - - - - - - - - - - - - - #

class aurora_parameters:
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

		self.fri_parameters = fri_parameters(\
			self.variables,					# Variable number (n)					\
			self.constraints,				# Constraints number (m)				\
			self.field_dim,					# Field dimension (fd)					\
			self.security_parameter + 3,	# Interactive Soundness Error (isp)		\
			self.security_parameter + 1, 	# Query Soundness Error (qsp)			\
			self.security_parameter*2, 		# Hash size (h)							\
			snd_type = snd_type, other_oracles = self.oracles_number())

	def __str__(self):
		return str(self.fri_parameters)

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

		fri_cost = self.fri_parameters.optimal_cost()
		extra_cost = self.extra_communication()
		return fri_cost + extra_cost
