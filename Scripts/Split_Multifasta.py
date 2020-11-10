#!/usr/bin/python
"""
Author: Anja Spang
takes multifasta file and splits it into several files with a desired amount of sequences
"""

import sys
import os
import argparse
from Bio import SeqIO
from collections import Counter

parser = argparse.ArgumentParser(description='extracts specific amount of sequences from a multifasta-file')
parser.add_argument('-m','--input1', metavar='FASTA', required = True, help='multifasta file from which sequences will be removed')
parser.add_argument('-n','--input2', metavar='int', required = True, help='integer: number of sequences per file')
args = parser.parse_args()


if args.input1 is None and not os.path.exists(args.input):
    print("input1 file missing!")
    sys.exit(-1)

if args.input2 is None and not os.path.exists(args.input):
    print("input2 missing!")
    sys.exit(-1)


multifasta = args.input1
seqnb = int(args.input2) - 1
split_count = 0
total_count = 0
filenb = 1


for seq_record in SeqIO.parse("%s" % (multifasta),"fasta"):
      header = str(seq_record.description).strip()
      seq = str(seq_record.seq).strip()
       
      if total_count == 0:
	  outfile = open("File%i.faa" % (total_count), "w")
     
      if split_count < seqnb:
	  split_count = split_count +1
	  outfile.write(">%s\n%s\n" % (header, seq))
      
      elif split_count == seqnb:
	  split_count = 0
	  outfile.close()
	  outfile = open("File%i.faa" % (total_count), "w")
	  outfile.write(">%s\n%s\n" % (header, seq))
	  filenb = filenb + 1
     
      total_count = total_count + 1

outfile.close()



print("Total amount of sequences: %s split into %s files" %(total_count, filenb))
		

	      
