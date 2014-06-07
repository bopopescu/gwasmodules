#!/usr/bin/env python

import os, sys
inputFname = os.path.expanduser("~/RefGenomes/dustPlus10_M1-22XY.bed.gz")
inputFname = os.path.expanduser("~/script/varcmp/scripts/LCR-hs37d5.bed.gz")
inputFname = os.path.expanduser("~/RefGenomes/dust_M1-22XY.bed.gz")
inputFname = os.path.expanduser("/illumina/scratch/CompetitiveAnalysis/CAG/Data/AnnotDB/Repeats/SegDups/genomicSuperDups_hg19.bed")

inputFname = os.path.expanduser("~/RefGenomes/dustPlus10_M1-22XY.overlap.genomicSuperDups_hg19.merged.bed")
inputFname=sys.argv[1]
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule import utils
from pymodule import MatrixFile
reader = MatrixFile(inputFname=inputFname, openMode='r', delimiter='\t')
span=0

for row in reader:
    if row[0][0]=='#':
        continue
    subSpan = int(row[2])-int(row[1]) + 1
    span += subSpan

print("span is %s \n"%(span))

