#!/usr/bin/python
#Anja, August 2017
#for python3.5
##########################################################################################
#Parse Blast results from diamond to keep only top blast hit, excluding self hit:
"""
depending on your blast output, something like this may work to select the columns of interest as shown in the example below:
awk -F'\t' -v OFS="\t" '{ print $1, $6, $7, $12, $13, $15, $16 }'
"""
#Has 7 columns separated by tabs (qseqid, salltitles, slen, evalue, bitscore, pident)
""""
accession 	s-description	s-length	E-value	Bitscore	identity taxID
CAP1_GB_A_00010 WP_013644295.1 diphthine synthase [Methanobacterium lacus]<>ADZ08944.1 diphthine synthase [Methanobacterium lacus]      263     2.4e-18 100.9   46.2    877455
CAP1_GB_A_00010 KYH41826.1 diphthine synthase [Candidatus Bathyarchaeota archaeon B26-2]        268     9.2e-18 99.0    42.0    1779371
CAP1_GB_A_00010 KYH42361.1 diphthine synthase [Candidatus Bathyarchaeota archaeon B26-1]        268     2.0e-17 97.8    42.0    1779370
CAP1_GB_A_00010 WP_071907154.1 diphthine synthase [Methanobacterium congolense]<>SCG86054.1 Diphthine synthase [Methanobacterium congolense]    263     5.0e-16 93.2    41.0    118062
CAP1_GB_A_00010 WP_069583342.1 diphthine synthase [Methanobacterium sp. A39]<>OEC87217.1 diphthine synthase [Methanobacterium sp. A39]  262     6.6e-16 92.8    38.5    1860100
"""


#Import modules
##########################################################################################

import sys
import os
import re
import argparse


#ARGPARSE; Open and read files
##########################################################################################

parser = argparse.ArgumentParser(description='Create file with one gene per line listing top hits')
parser.add_argument('-i','--input', metavar='TXT', required = True, help='Protein accessions of query sequences')
parser.add_argument('-d','--input2', metavar='tsv', required = True, help='Blast results using diamond')
parser.add_argument('-t','--input3', metavar='txt', required = False, help='if applicable: Taxa to be excluded, one taxon per line')
parser.add_argument('-o','--output', metavar='tsv', help='top blast hits', required = True)
args = parser.parse_args()

if args.input is None and not os.path.exists(args.input):
	print("input file missing!")
	sys.exit(-1)

if args.input2 is None and not os.path.exists(args.input2):
	print("input2 file missing!")
	sys.exit(-1)

if args.input3 is None:
	taxa_list = []
else:
	taxa_list =  open(args.input3).read().splitlines()


geneIDs = open(args.input).read().splitlines()
infile = open(args.input2).readlines()
outfile= open("%s" % args.output,'w')



#Functions to read in data
##########################################################################################


#Read Blast results into dictionary
def	Blast_result_all(diamond_blast):
	Blast_result_dic = {}

	for line in diamond_blast:
		entry = line.split('\t', 1) #splits only the first tab
		if entry[0] not in entry[1]:
			if entry[0] not in Blast_result_dic:
				Blast_result_dic[entry[0]] = entry[1].strip()
			else:
				pass

	return Blast_result_dic




#Read Blast results into dictionary excluding hits to selected taxa
def	Blast_result(diamond_blast, taxa_to_exclude):
	Blast_result_dic = {}

	for line in diamond_blast:
		entry = line.split('\t', 1)
		if entry[0] not in entry[1]:
			if entry[0] not in Blast_result_dic:
				if not any(tax in entry[1] for tax in taxa_to_exclude):
					Blast_result_dic[entry[0]] = entry[1].strip()
			else:
				pass

	return Blast_result_dic



#Parse and write final table:
##########################################################################################
"""
accession 	s-description-all	s-length-all	E-value	Bitscore-all	ID 		TaxID (s-description	s-length	E-value	Bitscore Identity	TaxID)
last 5 columns optional and only printed in case of input3 is true (i.e. if file with taxa to exclude is given)
"""



def Parse_final_table1(genes, Blast_all_dic):
	for geneID in genes:
		if geneID in Blast_all_dic:
			Blast_all = Blast_all_dic[geneID]
		else:
			Blast_all = "-\t-\t-\t-\t-\t-"

		outfile.write("%s\t%s\n" % (geneID, Blast_all))


def Parse_final_table2(genes, Blast_all_dic, Blast_sel_dic):
	for geneID in genes:
		if geneID in Blast_all_dic:
			Blast_all = Blast_all_dic[geneID]
		else:
			Blast_all = "-\t-\t-\t-\t-\t-"

		if geneID in Blast_sel_dic:
			Blast_sel = Blast_sel_dic[geneID]
		else:
			Blast_sel = "-\t-\t-\t-\t-\t-"

		outfile.write("%s\t%s\t%s\n" % (geneID, Blast_all, Blast_sel))




#Writing final table
##########################################################################################

def write_table1():
		outfile.write("accession\tTop-hit\tslenght\tE-value\tScore\tID%\tTaxID\n")
		Parse_final_table1(geneIDs, Blast_result_all(infile))

def write_table2():
		outfile.write("accession\tTop-hit\tslenght\tE-value\tScore\tID%\tTaxID\tTop-hit_wo\tslenght\tE-value\tScore\tID%\tTaxID\n")
		Parse_final_table2(geneIDs, Blast_result_all(infile), Blast_result(infile, taxa_list))

if __name__ == "__main__":
	if not taxa_list:
		write_table1()
	else:
		write_table2()



#Close files
##########################################################################################
outfile.close()
