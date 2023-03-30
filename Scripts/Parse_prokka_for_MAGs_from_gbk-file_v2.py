#!/usr/bin/python
# written by Anja Spang and adjusted by Nina Dombrowski
# The goal of this script is to extract protein features from genbank files

from Bio import SeqIO
import sys
import gzip
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re


parser = argparse.ArgumentParser(description='extract protein features from many separate gzipped genbank files')
parser.add_argument('-i','--input', metavar='TXT', required = True, help='text file giving the names of all zipped gbk files per lines from which to extract the protein sequences')
parser.add_argument('-t','--output1', metavar='tsv', help='tab file with MAG name, contig name, contig length in bp, protein accession, and corresponding annotations', required = True)
args = parser.parse_args()


if args.input is None and not os.path.exists(args.input):
  print("input file missing!")
  sys.exit(-1)


infile_list=open(args.input).readlines()
output_handle1 = open(args.output1,'w')


for file in infile_list:
	for seq_record in SeqIO.parse(gzip.open(file.strip()),"genbank"):
		#MAG_contig_len = seq_record.annotations['LOCUS']
		#contig_len= seq_feature.source.split('..')[1]
		#print(contig_len)
		contig = seq_record.id
		contig_len = len(seq_record.seq)		

		for seq_feature in seq_record.features:
			
			if seq_feature.type=="CDS" and 'pseudo' not in seq_feature.qualifiers:

				MAG= seq_feature.qualifiers['locus_tag'][0].rsplit('_',1)[0] if "locus_tag" in seq_feature.qualifiers else "-"

				locus_tag= seq_feature.qualifiers['locus_tag'][0] if "locus_tag" in seq_feature.qualifiers else "-"
				
				product= seq_feature.qualifiers['product'][0] if "product" in seq_feature.qualifiers else "-"
				
				start = seq_feature.location.start
				end = seq_feature.location.end
				strand = seq_feature.strand
				
				prokka_annotation = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (MAG, contig, contig_len, locus_tag, start, end, strand,  product)
				

				output_handle1.write("%s\n" % (prokka_annotation))





output_handle1.close()
