from .constants import *
from .utilities import *

import schemes as sch

# - - - - - - - - - - - - - - - #
# 	Testing functionalities 	#
# - - - - - - - - - - - - - - - #

def _check_file(filename):
	if os.path.isfile(filename):
		user_input = str_input(str_bold("Warning: ") + "\"{}\" already exists. Overwrite? (y/n) ".format(filename))
		
		if user_input.lower() not in ["y", "yes"]:
			print("Aborting.")
			sys.exit(0)

def _parameters_class(scheme):
	# Helper function for comparison_to_csv
	#  return the parameter class requested from the test
	#  
	#  scheme : "ligero" or "aurora"

	if scheme == "ligero":
		return sch.ligero_parameters
	elif scheme == "aurora":
		return sch.aurora_parameters
	else:
		# debug
		assert False, "Requested unkown parameters class"

def comparison_to_csv(scheme, filename, dim_min, dim_max, fd, sp, rmfe_iter, snd_type):
	# Make a comparison between plain scheme over F2 and our optimised version
	#
	#  scheme 		: "aurora" or "ligero" - choose which scheme to test
	#  filename 	: name of the file in which final data is stored
	#  dim_min 		: input size in the test start from 2**dim_min
	#  dim_max 		: input size in the test end with 2**dim_max (inclusive)
	#  fd 			: Aurora's field dimension
	#  rmfe_iter 	: A list of tuples (k,m) of rmfe parameters used in the optimised version
	#  sp 			: security parameter
	#  snd_type 	: soundness type

	filename = filename + ".csv"
	_check_file(filename)

	print_message("Storing results in {}\n".format(filename))
		
	# out : list of vectors (= columns) later copied in the csv
	out = []
	vec_constraints_dim = [n for n in range(dim_min, dim_max + 1)]
	out.append(["constraints"] + [2**n for n in vec_constraints_dim])

	# choose the parameter class
	parameters_class = _parameters_class(scheme)

	try:
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
		print_separator(verbose = True)
		out.append(["standard"] + vec_results)

		for rmfe in rmfe_iter:
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
			print_verbose_message("End of optimised {:s} tests, rmfe = {}".format(scheme, rmfe))
			print_separator(verbose = True)
			out.append(["optimised rmfe {:3d} {:3d}".format(*rmfe)] + vec_results)

	finally:
		# write to csv file
		with open(filename, "w") as file:
			vec_to_csv(file, out)

def general_test(scheme, dim_min, dim_max, fd, sp, rmfe_iter, snd_type):
	# General test suite. Return all the parameters found during the optimisations process
	#
	# scheme 	: "aurora" or "ligero" or other implemented schemes (see -help)
	# dim_min 	: input size in the test start from 2**dim_min
	# dim_max 	: input size in the test end with 2**dim_max (inclusive)
	# fd 		: Aurora's field dimension
	# rmfe 		: A tuple (k,m) of rmfe parameters used in the opt
	# sp 		: security parameter
	# snd_type 	: soundness type	parameters_class = _parameters_class(scheme)
	
	parameters_class = _parameters_class(scheme)

	for dim in range(dim_min, dim_max + 1):
		print_message("Testing standard {:s}, n = {:2d}".format(scheme, dim))

		# initialize the parameters of the selected scheme
		parameters = parameters_class(2**dim, 2**dim, fd, sp, snd_type, "standard")

		# find optimal parameters and get the nimum cost associated
		cost = parameters.optimal_cost()

		# print the optimal parameters and the cost found
		print_message(str(parameters))
		print_cost("Final cost", cost)
		print_separator(verbose = True)

	for rmfe in rmfe_iter:
		for dim in range(dim_min, dim_max + 1):
			print_message("Testing optimised {:s}, n = {:2d}, rmfe = {}".format(scheme, dim, rmfe))

			# initialize the parameters of the selected scheme
			parameters = parameters_class(2**dim, 2**dim, fd, sp, snd_type, "optimised", rmfe = rmfe)

			# find optimal parameters and get the nimum cost associated
			cost = parameters.optimal_cost()

			# print the optimal parameters and the cost found
			print_message(str(parameters))
			print_cost("Final cost", cost)
			print_separator(verbose = True)

	print_message("End of test")

