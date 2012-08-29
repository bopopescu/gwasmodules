#!/usr/bin/env python
""""

Examples:
	%s
	
	%s -u yh -p secret -d vervetdb -k taxonomy

Description:

2012.8.28 test class for TaxonomyDB.py
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import unittest, os, sys, getopt, csv
from pymodule.TaxonomyDB import TaxonomyDB

class TestTaxonomyDB(unittest.TestCase, TaxonomyDB):
	__doc__ = __doc__
	option_default_dict = TaxonomyDB.option_default_dict.copy()
	def __init__(self, testname=None, **keywords):
		"""
		2012.8.28
		"""
		TaxonomyDB.__init__(self, **keywords)
		
		super(TestTaxonomyDB, self).__init__(testname)
		#self.hs = helpspot.HelpSpot(path, user, pword)
		
	def setUp(self):
		print
	
	def test_returnScientificNameGivenTaxID(self):
		self.setup()
		taxID = 60711
		
		#1st just db query
		scientificName = self.getScientificNameGivenTaxID(taxID)
		self.assertEqual(scientificName, "Chlorocebus sabaeus")
		
		#2nd cached version
		scientificName = self.returnScientificNameGivenTaxID(taxID)
		sys.stderr.write("Scientific name for %s is %s.\n"%(taxID, scientificName))
		self.assertEqual(scientificName, "Chlorocebus sabaeus")

	
	def test_returnTaxIDGivenScientificName(self):
		self.setup()
		scientificName = "Chlorocebus sabaeus"
		#1st just db query
		taxID = self.getTaxIDGivenScientificName(scientificName)
		self.assertEqual(taxID, 60711)
		
		#2nd cached version
		taxID = self.returnTaxIDGivenScientificName(scientificName)
		sys.stderr.write("Taxonomy ID for %s is %s.\n"%(scientificName, taxID))
		
		self.assertEqual(taxID, 60711)

		
	
if __name__ == '__main__':
	import sys, os, math
	
	TestCaseDict = {1:TestTaxonomyDB,}
	type = 1
	
	from pymodule import ProcessOptions
	main_class = TestCaseDict[type]
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	suite = unittest.TestSuite()
	suite.addTest(TestCaseDict[type](testname='test_returnScientificNameGivenTaxID', **po.long_option2value))
	suite.addTest(TestCaseDict[type](testname='test_returnTaxIDGivenScientificName', **po.long_option2value))
	#suite.addTest(unittest.makeSuite(TestCaseDict[type]))
	unittest.TextTestRunner(verbosity=2).run(suite)
	#unittest.TextTestRunner().run(suite)
	"""

	instance = main_class(**po.long_option2value)
	instance.setup()
	instance.run()
	"""
