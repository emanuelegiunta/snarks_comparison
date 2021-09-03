import math

from .constants import *


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
			sample = set(rnd.randint(0, int(leaf_num), int(delta)))
			leaves_set = leaves_set.union(sample)

	out = len(leaves_set)
	size_new = out

	k = int(math.ceil(np.log2(leaf_num)))
	f = lambda i, n : n % (1 << (k - i))
	
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
	# Evaluate the communication cost of BSC, i.e. the cost of opening
	#  oracles, made of three terms: the size of acepting paths in the 
	#  associated Merkel Tree and the actual message sent, an element
	#  of the alphabet. Finally (rounds + 1) seeds - obtained as the hash
	#  of previous messages - are required.
	#
	# hash_size 		: bit length of the hash output used (128/256)
	# vec_alphabet_size : bit size of the alphabet used for the oracles (i.e. size of the field)
	# rounds 			: number or rounds
	# vec_oracle_length : length of the oracles sent
	# vec_query 		: number of queries to the i-th oracle

	# adapt single values for vec_oracle_length and vec_query:
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
		# compute the number of mkt hashes required from Monte-Carlo evaluation
		mkt_opening = hash_size*sum(mkt_avg_cost(vec_oracle_length[i], vec_query[i]) for i in range(len(vec_query)))
		mkt_opening = int(math.ceil(mkt_opening))
	else:
		# compute the number of mkt hashes required from an estimate of the expected value, faster
		mkt_opening = hash_size*sum(mkt_estimate_cost(vec_oracle_length[i], vec_query[i]) for i in range(len(vec_query)))

	# cost of the opened items
	mkt_opened_elements = sum(vec_alphabet_size[i] * vec_query[i] for i in range(len(vec_query))) 
	
	# cost of the seed for verifiers replies
	mkt_round_hash = hash_size*(rounds + 1) 

	out = mkt_opening + mkt_opened_elements + hash_size*(rounds + 1) 
	return out