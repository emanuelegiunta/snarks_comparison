import numpy as np
import numpy.random as rnd
import math
import matplotlib.pyplot as plt
import sys 
import os.path

from .constants import *

class _cst:
	def __init__(self):
		pass

cst = _cst()

cst.VERBOSE_FLAG = VERBOSE_FLAG
cst.VERY_VERBOSE_FLAG = VERY_VERBOSE_FLAG

class CompatibilityError(Exception):
	pass

class InputError(Exception):
	pass


# - - - - - - - - - - - - - - - - - - - #
#	basic math and generic functions	#
# - - - - - - - - - - - - - - - - - - - #

def comb(n,m):
	'''binomial coefficient of n, m

	i.e. n! / ((n - m)! * m!)
	'''
	out = 1
	for i in range(m):
		out *= (n - i)/float(m - i)
	return out

def combr(n, m, k):
	'''ratio of comb(n,k)/comb(m,k)

	Note that it is faster and safer to run combr(m,k) instead of evaluating 
	both numerator and denominator, as they grow exponentially in k
	'''
	out = 1
	for i in range(int(k)):
		out *= (n - i)/float(m - i)
	return out

def xor(a, b):
	''' xor operator

	if a and b are functions' outputs, those functions are evaluated exactly
	once
	'''

	return (a or b) and not (a and b)

def _ceil_decorator(fun):
	'''ceil decorator - redefine ceil as math.ceil if the latter follows the
	type signature

	ceil(n: Union[int, float]) -> int
	'''
	if isinstance(math.ceil(.0), int):
		return math.ceil 
	else:
		return fun

@_ceil_decorator
def ceil(n):
	return int(math.ceil(n))

def ceilpow2(n):
	''' return the smaller integer greater than n that is a power of 2
	
	Note: a quick implementation instead of an efficient one is currently used
	'''

	return 2**int(math.ceil(np.log2(n)))

def binsearch(f, x_min, upper_bound = None):
	'''return the minimum x such that f(x) is True
	f 		: a boolean valued function. For the bisection to make sense, f
			: should be "monotone" in {False, True}
	x_min 	: minimum value of x. we assume x_min > 0
	'''

	if upper_bound is None:
		upper_bound = float('inf')

	# Find an upper bound
	x_max = 2*x_min

	while not f(x_max):
		# abort if no solution lies in the given range
		if x_max > upper_bound:
			return float('inf')

		x_max *= 2

	# next we run the binary search to minimize f
	while(x_max - x_min > 1):
		if f((x_min + x_max)//2) :
			x_max = (x_min + x_max)//2
		else:
			x_min = (x_min + x_max)//2

	return x_max

# - - - - - - - - - - - - - - - - - - - - - #
#	String formatting and message handling	#
# - - - - - - - - - - - - - - - - - - - - - #

def print_cost(name, bit_cost):
	if cst.VERBOSE_FLAG:
		print((name + ":\t%10d bit%10d B%10d KB") % (bit_cost, bit_cost/8, bit_cost/8000) )

def print_verbose_cost(name, bit_cost):
	if cst.VERY_VERBOSE_FLAG:
		print_cost(name, bit_cost)

def print_message(string):
	if cst.VERBOSE_FLAG:
		print(string)

def print_verbose_message(string):
	'''print a message if the VERY_VERBOSE_FLAG is true
	string	: input message
	'''

	if cst.VERY_VERBOSE_FLAG and cst.VERBOSE_FLAG:
		print(string)

def print_separator(length = 36, indentation_length = 2, verbose = False, very_verbose = False):
	'''print a line to separate data
	length 				: number of caracter printed
	indentation_length 	: number of spaces before the line
	very_verbose 		: if True check the cst.VERY_VERBOSE_FLAG
	'''

	if (very_verbose and cst.VERY_VERBOSE_FLAG) or (verbose 
		and cst.VERBOSE_FLAG):

		# print line
		print("\n" + " "*indentation_length + "="*length + "\n")

def vec_to_csv(file, vec):
	'''Given a list of list, print each list as a csv column in file
	file 	: file object, already opened
	vec 	: list of list of the same length
	'''

	# debug
	assert len( {len(col) for col in vec} ) <= 1, "Writing columns of different length in csv"

	#werite a vector into an open csv file
	for i in range(len(vec[0])):
		file.write(",".join([repr(a[i]) for a in vec]) + "\n")

def str_bold(s):
	'''print bold text, need the terminal to support ANSI formatting
	'''

	return "\033[1m{:s}\033[0m".format(s)

def str_underline(s):
	'''print underlined text, need the terminal to support ANSI formatting
	'''

	return "\033[4m{:s}\033[0m".format(s)

def str_synopsis():
	'''Return the signature of the main programm upon receiving ill formed input
	'''

	out  = ""
	out += str_bold("SYNOPSIS: ") + "{}\n".format(sys.argv[0])
	out += "\t[-l] [-tl] [-a] [-ta] [-sp {arg}] [-ld {arg}] [-hd {arg}]\n"
	out += "\t[-h] [-p] [-r {arg}] [-fd {arg}] [-f {arg}]\n"
	out += "\t[-D {arg}] [-v] [-vv] [-help]\n"
	out += "\n"
	
	out = out.format(arg = str_underline("1 argument"), args = str_underline("2 arguments"))
	return out

def str_help():
	'''Return the message printed with the -help flag in any position
	'''

	out  = ""
	out += str_bold("Parameter optimiser for Aurora [EC:BCRSVW19], "
		"Ligero [CCS:AHIV17], Ligero++ [CCS:BFHVXZ20] with boolean R1CS\n\n")
	out += str_synopsis()

	out += str_bold("COMMAND LINE DESCRIPTION:\n")
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-l, --ligero"))
	out += " Runs a comparison between standard and optimised Ligero\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-a, --aurora"))
	out += " Runs a comparison between standard and optimised Aurora\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-lpp, --ligeropp"))
	out += " Runs a comparison between standard and optimised Ligero++\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-tl, --test_ligero"))
	out += " Runs a test for Ligero, printing used parameters on screen\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-ta, --test_aurora"))
	out += " Runs a test for Ligero++, printing used parameters on screen\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-tlpp, --test_ligeropp"))
	out += " Runs a test for ligero, printing used parameters on screen\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-sp, --security_parameter"), str_underline("(1 argument)"))
	out += " Set the security parameter\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-ld, --lowest_dimension"), str_underline("(1 argument)"))
	out += " Set the minimum size of tested constraints, inclusive\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-hd, --highest_dimension"), str_underline("(1 argument)"))
	out += " Set the maximum size of tested constraints, inclusive\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-h, --heuristic"))
	out += " Uses optimistic soundness bounds\n"
	out += "\n"

	out += "  {:s}:\n".format(str_bold("-p, --proven"))
	out += " Uses proven soundness bounds, currently the default option\n"
	out += "\n"

	out += "  {:s} {:s}:\n".format(str_bold("-r, --rmfe"), str_underline("(1 argument)"))
	out += (" Set the reverse multiplication-friendly embedding (k,m) used by "
		"the optimised version. Input should countain at least one ore more "
		"couples of integers and has to be placed in quotes of double quotes "
		"and separated by commas, spaces or semicolons (or a combination of "
		"them). Brakets are ignored. Below examples of valid input\n"
		"\t\"48 192\" parsed as [(48, 192)]\n"
		"\t\'[3,5];[2,3]\' parsed as [(3,5), (2,3)]\n"
		"\t\"[(42; 135), [2 3]]\" parsed as [(42, 135), (2, 3)]\n")
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

	return out

def str_input(s = ""):
	'''Prompt the user with s and reads the input on the stdin
	
	Compatible with both python 2 and 3
	'''
	
	if sys.version_info.major == 2:
		return raw_input(s)
	else:
		return input(s)
