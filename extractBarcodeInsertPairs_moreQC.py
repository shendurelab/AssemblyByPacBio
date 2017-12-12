#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *01.03.2017
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 
import pysam

def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op == 4: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def parseMD(MD):
  MDfields = []
  value = ""
  chars = ""
  for elem in MD:
    if elem.isdigit() and chars == "":
      value+=elem
    elif not elem.isdigit() and value == "":
      chars+=elem
    elif elem.isdigit() and chars != "":
      MDfields.append(chars)
      chars = ""
      value = elem
    elif not elem.isdigit() and value != "":
      MDfields.append(int(value))
      value = ""
      chars = elem
  if value != "": MDfields.append(int(value))
  if chars != "": MDfields.append(chars)
  return MDfields

def parseMDwithCigar(MD,cigar,seq):
  MDfields = parseMD(MD)
  res = []
  posRef = 0
  posRead = 0
  ind = 0 # Pointer in MDfields
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    #print op,ind
    if (op in [0,7,8]) and ind % 2 == 1:
      while (count > 0):
        for base in MDfields[ind]:
          res.append((posRef,posRead,(seq[posRead],base)))
          posRef += 1
          posRead += 1
          count-= 1
        ind += 1 # Next posRef should be a number
        if (count > 0):
          if MDfields[ind] > count:
            posRef += count
            posRead += count
            MDfields[ind] -= count
            count = 0
          else:
            count -= MDfields[ind]
            posRef += MDfields[ind]
            posRead += MDfields[ind]
            ind += 1
      if (ind < len(MDfields)) and (MDfields[ind] == 0):
        ind += 1
    elif (op in [0,7,8]) and ind % 2 == 0:
      if MDfields[ind] > count:
        posRef += count
        posRead += count
        MDfields[ind] -= count
      else:
        count -= MDfields[ind]
        posRef += MDfields[ind]
        posRead += MDfields[ind]
        ind += 1 # Next posRef should be a base/string of bases
        while (count > 0):
          for base in MDfields[ind]:
            res.append((posRef,posRead,(seq[posRead],base)))
            posRef+= 1
            posRead+= 1
            count-= 1
          ind += 1 # Next posRef should be a number
          if (count > 0):
            if MDfields[ind] > count:
              posRef += count
              posRead += count
              MDfields[ind] -= count
              count = 0
            else:
              count -= MDfields[ind]
              posRef += MDfields[ind]
              posRead += MDfields[ind]
              ind += 1
      if (ind < len(MDfields)) and (MDfields[ind] == 0):
        ind += 1
    elif (op == 4): 
      posRead += count
    elif (op == 5):
      continue
    elif (op == 1):
      res.append((posRef,posRead,("INS",count,seq[posRead:posRead+count])))
      posRead += count
    elif (op == 2) and (ind % 2 == 0):
      posRef += count
      ind+=2
    elif (op == 2) and (ind % 2 == 1):
      res.append((posRef,posRead,("DEL",len(MDfields[ind][1:]),MDfields[ind][1:])))
      posRef += count
      ind += 1
    else:
      sys.stderr.write("Cigar case not implemented: (%d,%d) (%s) [%d (%s),%d]\n"%(op,count,str(cigar),ind,str(MDfields),posRef))
      #sys.exit()
      return []
  #print MD,cigar,MDfields,res,posRef
  return res

def extractRegionRefCoords(offset, alignmentDiffs, seq, qual, rstart, rend):
  rstart-=(offset-1)
  rend-=(offset-1)
  #start,end = 0,0
  start,end = rstart,rend
  for posRef,posRead,event in alignmentDiffs:
    if posRef <= rstart: start = posRead+(rstart-posRef)
    if posRef <= rend: end = posRead+(rend-posRef)
    #print start,end
    if (event[0] == "INS"): 
      posRead += event[1]+1
      posRef += 1
      if posRef <= rstart: start = posRead+(rstart-posRef)
      if posRef <= rend: end = posRead+(rend-posRef)
    elif (event[0] == "DEL"): 
      posRead += 1
      posRef += 1
      for i in range(event[1]):
        posRef += 1
        if posRef <= rstart: start = posRead+(rstart-posRef)
        if posRef <= rend: end = posRead+(rend-posRef)
    else: continue
    #print start,end
    
  if (start < 0) or (end > len(seq)): return None,None
  else: return seq[start:end],qual[start:end]
  

parser = OptionParser()
parser.add_option("-l","--length", dest="blength", help="Length of the barcodes (default 16)",default=16,type="int")
parser.add_option("-p","--position", dest="bpos", help="Position of the barcode (default 52)",default=52,type="int")
parser.add_option("--start", dest="start", help="Start reference position for extracted sequence (default 114)",default=114,type="int")
parser.add_option("--end", dest="end", help="End reference position for extracted sequence (default 1065)",default=1065,type="int")

parser.add_option("-b","--minBaseQ", dest="minBaseQ", help="Minimum base quality score (default 0)",default=0,type="int")
parser.add_option("-q","--minMapQ", dest="minMapQ", help="Filter alignments by mapQ (default 0)",default=0,type="int")
parser.add_option("-s","--allow-soft-clipped", dest="allowSoftClipped", help="Consider soft clipped reads (default Off)",default=False,action="store_true")
#parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default reference.fa)",default="reference.fa")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

#genome = pysam.Fastafile(options.reference)

count = 0
passedBasic, extracted, correctTagLength = 0,0,0
dup = 0
qcfail = 0
unmapped = 0
badqual = 0
softClipped = 0
#options.bpos = options.bpos-1
options.start = options.start-1
options.end = options.end-1

for bamfile in args:
  if options.verbose:
    print "Opening:",bamfile
  input_file = pysam.Samfile( bamfile, "rb" )
  BAMreferences = dict(enumerate(input_file.references))
  #referenceLengths = dict(zip(input_file.references,input_file.lengths))
  #if options.verbose:
    #print str(BAMreferences)[:200]
    #print str(referenceLengths)[:200]

  chrom = None
  for read in input_file:
    count += 1

    
    if read.qual == None or len(read.qual) != len(read.seq):
      read.qual="!"*len(read.seq)
      badqual += 1
    if read.is_duplicate:
    	 dup +=1
    	 continue
    	 
    if read.is_qcfail :
    	qcfail +=1
    	continue
    	
    if read.is_unmapped:
    	unmapped +=1
    	continue
    	
    if isSoftClipped(read.cigar) and not options.allowSoftClipped:
    	softClipped += 1
    	continue
    if options.minMapQ > read.mapq: continue
    if options.minBaseQ > min(map(lambda x: ord(x)-33,read.qual)): continue
    passedBasic += 1
    #if count > 25: sys.exit()
    
    cchrom = BAMreferences[read.tid]
    rstart,rend = 0,0
    if read.is_paired:
      if read.mate_is_unmapped: continue
      if read.rnext != read.tid: continue
      rstart = read.pos # 0-based
      rend = rstart+aln_length(read.cigar) # end excluded
    else:
      rstart = read.pos # 0-based
      rend = rstart+aln_length(read.cigar) # end excluded

    #sys.stderr.write(read.qname+"\n")
    for key,value in read.tags:
      if key == "MD":
        alignmentDiffs = parseMDwithCigar(value,read.cigar,read.seq)
        for posRef,posRead,event in alignmentDiffs:
          if (event[0] == "INS") and (event[1] == options.blength) and (options.bpos-5 <= posRef+rstart <= options.bpos+5): 
            correctTagLength += 1
            seq,qual = extractRegionRefCoords(read.pos, alignmentDiffs, read.seq, read.qual, options.start, options.end)
            if seq != None:
              sys.stdout.write(read.seq[(posRead-(posRef+rstart-options.bpos)):(posRead-(posRef+rstart-options.bpos)+event[1])]+"\t%s\t%s\n"%(seq,qual)) # read.qname,posRef,posRead,event,"\t"+event[2],"\t"+"%d,%d"%(posRef,rstart)
              extracted += 1

  input_file.close()

if options.verbose:
  sys.stderr.write("Read %d alignments, %d badqual,  %d dups, %d qcfail, %d unmapped, %d softClipped,  %d passed filters, %d insertion matched barcode length, %d reported back\n"%(count,badqual,dup,qcfail,unmapped,softClipped,passedBasic,correctTagLength,extracted))