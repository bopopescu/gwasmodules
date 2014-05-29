#!/usr/bin/env python
"""
2012.11.20 data structures for IO
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from HDF5MatrixFile import HDF5MatrixFile
from SNP import getGenomeWideResultFromHDF5MatrixFile
from Association import constructAssociationPeakRBDictFromHDF5File, getAssociationLandscapeDataFromHDF5File