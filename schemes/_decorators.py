def protocol_types(*argv):
	
	# Create the parametrised decorator
	def _protocol_types(cls):
		# Reject if "protocol_types" was defined within cls
		if hasattr(cls, "protocol_types"):
			raise RuntimeError("{} has already protocol_types attribute".format(
				cls.__name__))

		# Store the types name in "protocol_types"
		#  by default always have "standard"
		cls.protocol_types = ["standard"] + [s for s in argv]
		return cls

	# If there is only one argument, that is not a string, interpret it as a
	#  class
	if len(argv) == 1 and not isinstance(argv[0], str):
		# _protocol_types receives no argument in this case
		input_class = argv[0]
		argv = []
		return _protocol_types(input_class)

	else:
		return _protocol_types
