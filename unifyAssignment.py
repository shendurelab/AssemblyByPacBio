#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *01.03.2016
"""

import sys, os
import gzip
from optparse import OptionParser
from collections import defaultdict

parser = OptionParser()
#parser.add_option("-b","--barcodes", dest="barcodes", help="File with barcodes (default 'input.txt')",default="input.txt")
#parser.add_option("-m","--minQualityScore", dest="minQual", help="Minimum quality threshold (default 40)",default=40,type="int")
#parser.add_option("-o","--numberOfOutliers", dest="QualOutliers", help="Maxmimum number of bases below quality threshold (default 0)",default=0,type="int")
#parser.add_option("-c","--minBarcodeCount", dest="minCount", help="Minimum barcode count (default 5)",default=5,type="int")
(options, args) = parser.parse_args()

def mean(numbers): 
  if len(numbers) > 0:
    return sum(numbers)/float(len(numbers))
  else:
    return 0

barcodes = defaultdict(list)

for filename in args:
  if os.path.exists(filename):
    if filename.endswith(".gz"):
      infile = gzip.open(filename)
    else:
      infile = open(filename)
    for line in infile:
      fields = line.split()
      if len(fields) == 3:
        mQual = mean(map(lambda x:ord(x)-33,fields[2]))
        barcodes[fields[0]].append((fields[1],mQual))
    infile.close()

for barcode,observations in sorted(barcodes.iteritems()):
  freq = defaultdict(int)
  for seq,qual in observations:
    freq[seq]+=1
  res = map(lambda (x,y):(y,x),freq.iteritems())
  res.sort()
  if res[-1][0]/float(len(observations)) > 0.5:
    sys.stdout.write("%s\t%s\n"%(barcode,res[-1][1]))
  else:
    res = map(lambda (x,y):(y,x),observations)
    res.sort()
    sys.stdout.write("%s\t%s\n"%(barcode,res[-1][1]))