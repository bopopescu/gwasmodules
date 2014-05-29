import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
    sys.path.insert(0, os.path.expanduser('~/lib64/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
try:
	import psycopg2 as psycopg
except ImportError:
	try:
		import psycopg
	except ImportError:
		sys.stderr.write("Neither psycopg nor psycopg2 is installed.\n")
import sys, getopt, csv, re
#from codense.common import db_connect, dict_map
from pymodule.yhio.SNP import ab2number, number2ab, NA_set, nt2number, \
	number2nt, number2color, nt_number_matching_matrix, get_nt_number2diff_matrix_index,\
	RawSnpsData_ls2SNPData, transposeSNPData, SNPData2RawSnpsData_ls
from sets import Set
#import networkx as nx
#from pymodule import importNumericArray, ProcessOptions, SNPData, write_data_matrix, read_data, TwoSNPData
#num = importNumericArray()

from db import Stock_250kDB, StockDB, AtDB

#2012.11.18
from mapper.AbstractVariationMapper import AbstractVariationMapper
#2012.11.18
from pegasus.AbstractVariationWorkflow import AbstractVariationWorkflow
#2012.11.18 DefineAssociationLandscapePipeline has to be after AbstractVariationWorkflow
# because it's derived from AbstractVariationWorkflow.
from association_peak.DefineAssociationLandscapePipeline import DefineAssociationLandscapePipeline
