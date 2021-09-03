from util import *
from schemes import *


if __name__ == "__main__":
	# drop the program name
	argv = sys.argv[1:]

	# create a parser object that handles cmd's input
	parser = parser_input()

	# parse the input
	parser.parse(argv)

	# complete missing fields
	parser.complete()

	# debug - all the parser fields has to be filled
	assert parser.check_empty_entries(exclude = True, var = ['fd']) is None, \
		"Parser has empty variable {}".format(
		parser.check_empty_entries(exclude = True, var = ['fd']))
	# - - - - - end of debug - - - - - #

	# run the requested process
	parser.run()
