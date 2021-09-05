from .constants import *
from .utilities import *

from .tests import *

# - - - - - - - - - - - - - #
#	Input Parser and class 	#
# - - - - - - - - - - - - - #

def _parser_rmfe(s):
	# Returns a list if the given string is compliant
	#  (see -help to check current correct form)
	#  otherwise return None.

	# step 1: we strip brackets, commas and semicolons
	for c in "}])([{,;":
		s = s.replace(c, " ") # maybe we should replace this with a space

	# step 2: we collaps consecutive spaces into one
	while "  " in s:
		s = s.replace("  ", " ")

	# step 3: we remove trailing spaces
	if s[0] == " ":
		s = s[1:]
	if s[-1] == " ":
		s = s[:-1]

	# step 4: we split
	out = s.split(" ")

	# we exclude the cases in which nothing is passed like "(,)" or ""
	if len(out) == 0:
		return None
	
	# step 5: produce the tuples assuming the string has an even number of entries 
	#  and those entries are integers
	try:
		return [( int(out[2*i]), int(out[2*i + 1]) ) for i in range(len(out)/2)]
	except (ValueError, IndexError) as e:
		return None

class parser_input:
	def __init__(self):
		self.test_type = "compare"
		self.scheme = "aurora"
		self.sp = 128
		self.dim_min = 8
		self.dim_max = 20
		self.rmfe_iter = None 		#auto completed after parsing
		self.snd_type = "proven" 
		self.fd = None				#auto completed after parsing
		self.filename = None 		#auto completed after parsing
		self.dirname = ""

	def check_empty_entries(self, var = None, exclude = False):
		# check that all variables in var are not set to None.
		#  if exclude, then check all variables beside those
		#  in var
		# 
		# It returns the name of the empty variable. If no empty
		#  variable is found, it returns None
		if var is None:
			var = []

		for key, value in vars(self).items():
			if value == None and (key in var or exclude) and not (key in var and exclude):
				print(key)
				return key
		else:
			return None

	def parse(self, argv):
		#
		# -l 	--ligero 				Runs a comparison for Ligero
		# -a 	--aurora 				Runs a comparison for Aurora
		# -lpp	--ligero++				Runs a comparison for Ligero++
		#
		# -tl   --test_ligero			Runs a test for Ligero
		# -ta 	--test_aurora 			Runs a test for Aurora
		# -tlpp --test_ligero++			Runs a test for Ligero++
		#
		# -sp 	--security_parameter	Set the security parameter (1 arg)
		# -ld 	--lowest_dimension		Set the minimum dimension tested, inclusive (1 arg)
		# -hd 	--highest_dimension 	Set the maximum dimension tested, inclusive (1 arg)
		# -h 	--heuristic 			Use heuristic soundness bounds (0 arg)
		# -p 	--proven 				Use proven soudness bounds (0 arg)
		# -r 	--rmfe 					Set the RMFE for the optmised version (1 arg)
		# -fd 	--field_dimension 		Set the field dimension for the standard version (1 arg)
		# -f 	--file_name 			Choose the file name (1 arg)
		# -D 	--directory				Absolute or relative path, added to filename as a prefix (1 arg)
		# -v 	--verbose 				Print some information on the screen
		# -vv 	--very_verbose 			Print more informations on the screen

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
					self.test_type = "compare"
					self.scheme = "ligero"
				
				elif key in ["-a", "--aurora"]:
					self.test_type = "compare"
					self.scheme = "aurora"

				elif key in ["-lpp", "--ligero++"]:
					self.test_type = "compare"
					self.scheme = "ligero++"

				elif key in ["-tl", "--test_ligero"]:
					self.test_type = "test"
					self.scheme = "ligero"
					cst.VERBOSE_FLAG = True

				elif key in ["-ta", "--test_aurora"]:
					self.test_type = "test"
					self.scheme = "aurora"
					cst.VERBOSE_FLAG = True

				elif key in ["-tlpp", "--test_ligero++"]:
					self.test_type = "test"
					self.scheme = "ligero++"
					cst.VERBOSE_FLAG = True

				elif key in ["-sp", "--security_parameter"]:
					self.sp = int(argv[i+1])
					i += 1

				elif key in ["-ld", "--lowest_dimension"]:
					self.dim_min = int(argv[i+1])
					i += 1

				elif key in ["-hd", "--highest_dimension"]:
					self.dim_max = int(argv[i+1])
					i += 1

				elif key in ["-h", "--heuristic"]:
					self.snd_type = "heuristic"

				elif key in ["-p", "--proven"]:
					self.snd_type = "proven"

				elif key in ["-r", "--rmfe"]:
					self.rmfe_iter = _parser_rmfe( argv[i+1] )
					if self.rmfe_iter == None:
						raise InputError
					i += 1

				elif key in ["-fd", "--field_dimension"]:
					self.fd = int(argv[i+1])
					i += 1

				elif key in ["-f", "--file_name"]:
					self.filename = str(argv[i+1])
					i += 1

				elif key in ["-D", "--directory"]:
					self.dirname = str(argv[i+1])
					i += 1

				elif key in ["-v", "--verbose"]:
					cst.VERBOSE_FLAG = True

				elif key in ["-vv", "--very_verbose"]:
					cst.VERY_VERBOSE_FLAG = True
					cst.VERBOSE_FLAG = True

				else:
					# If there is no flag where a flag was expected
					#  raise an InputError 
					raise InputError

				i += 1
		except (IndexError, ValueError, InputError):
			print(str_synopsis())
			sys.exit(0)

	def complete(self):

		# set rmfe
		if self.rmfe_iter == None:
			if self.scheme == "aurora":
				self.rmfe_iter = [(48, 198)]
			elif self.scheme == "ligero":
				self.rmfe_iter = [(48, 160)]
			elif self.scheme == "ligero++":
				self.rmfe_iter = [(48, 160)]
			else:
				raise ValueError("unknown scheme requested")

		# set fd for those scheme that needs it
		if self.fd == None:
			if self.scheme == "aurora":
				self.fd = 192
			elif self.scheme == "ligero":
				pass
			elif self.scheme == "ligero++":
				self.fd = 128
			else:
				raise ValueError("unknown scheme requested")

		# set filename
		if self.filename == None:
			if len(self.rmfe_iter) == 1:
				self.filename = "{:s}_{:s}_sp_{:d}_rmfe_{:d}-{:d}".\
					format(self.scheme, self.snd_type, self.sp, *(self.rmfe_iter[0]))
			else:
				self.filename = "{:s}_{:s}_sp_{:d}_multi-rmfe".\
					format(self.scheme, self.snd_type, self.sp)

		# add the directory	
		self.filename = self.dirname + self.filename

	def run(self):
		# test execution
		if self.test_type == "compare":
			comparison_to_csv(self.scheme, self.filename, self.dim_min, self.dim_max, self.fd, self.sp, \
				self.rmfe_iter, self.snd_type)

		elif self.test_type == "test":
			general_test(self.scheme, self.dim_min, self.dim_max, self.fd, self.sp, self.rmfe_iter, self.snd_type)
