import numpy as np
import numpy.random as rnd
import itertools as it

# constants used by BSC
MKT_ITER = 1000					# recommended 1000	
MKT_ESTIMATE_ITER = 20			# recommended 20
MTK_FAST_FLAG = True 			# recommended True

# constants used by ligero.
LGR_MAX_SAMPLE = 150			# recommended 150
LGR_MAX_DOMAIN_DIM = 6			# recommended 6

# constants used by FRI
FRI_ITER = 1000					# recommended 1000
FRI_MAX_DOMAIN_DIM = 6			# recommended 6
FRI_MAX_LOCALIZATION_NUM = 3	# recommended 4

# constants used by Ligero++
LPP_MAX_DOMAIN_DIM = 6			# STILL TESTING
LPP_MAX_LOG_L = 10				# STILL TESTING