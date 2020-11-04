#!/usr/bin/python
#Anja, February 2017
#for python3.5
##########################################################################################
#Parse IPRscan result yielding a tab-delimited file with one gene per line and following columns separated  by tabs:
#GeneID	IPRdomain	IPRdescription	PFAMdomain	PFAMdescription	KEGGresults

#infile looks as following:
#accession in 1. column, PFAM in 5. column, PFAM description in 6. column, IPRdomain in 12th column, IPRdescription in 13th column, GO number in 14th column, EC number in 15th column.


"""Example
OLS14062        c04583dc6ecdecb405635bc83a3eb09a        159     Gene3D  G3DSA:2.40.50.140               8       155     8.3E-33 T       11-02-2017
OLS14062        c04583dc6ecdecb405635bc83a3eb09a        159     Pfam    PF01176 Translation initiation factor 1A / IF-1 29      103     5.5E-20 T       11-02-2017      IPR006196       RNA-binding domain, S1, IF1 type        GO:0003723|GO:0003743|GO:0006413
OLS14062        c04583dc6ecdecb405635bc83a3eb09a        159     SUPERFAMILY     SSF50249                12      150     1.8E-26 T       11-02-2017      IPR012340       Nucle

"""




##########################################################################################
#Import modules
##########################################################################################

from Bio import SeqIO
import sys
import gzip
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import fileinput



##########################################################################################
#ARGPARSE
##########################################################################################

parser = argparse.ArgumentParser(description='Create file with one gene per line listing all corresponding IPR domains, pfam domains and EC-numbers, (information extracted from IPR-result file')
parser.add_argument('-i','--input', metavar='TXT', required = True, help='Result file of IPR scan, tab formated')
parser.add_argument('-s','--input2', metavar='FASTA', required = True, help='Fasta_file with all proteins, fasta')
parser.add_argument('-o','--output', metavar='TXT', help='proteins with corresponding domain information', required = True)
args = parser.parse_args()

if args.input is None and not os.path.exists(args.input):
  print("input file missing!")
  sys.exit(-1)

if args.input2 is None and not os.path.exists(args.input):
    print("input2 file missing!")
    sys.exit(-1)



##########################################################################################
#Open and read files
##########################################################################################

infile = open(args.input).readlines()
#print(infile)
num_lines = sum(1 for line in open(args.input).readlines())
#print(num_lines)
infile2 = args.input2
outfile=open("%s" % args.output,'w')



##########################################################################################
#Functions
##########################################################################################

def make_list_geneID(fastafile):
	geneIDs = []

	for seq_record in SeqIO.parse("%s" % (fastafile),"fasta"):
		geneIDs.append(seq_record.id)

	return geneIDs



#Read genes with corresponding ipr domains into dictionary
##########################################################################################

def IPRdomain_dictionary(IPRscan_out):
	IPRdomains = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'IPR' in entry[11]:
				if entry[0] not in IPRdomains:
					IPRdomains[entry[0]] = []
					IPRdomains[entry[0]].append(entry[11])
				else:
					if entry[11] not in IPRdomains[entry[0]]:
						IPRdomains[entry[0]].append(entry[11])
		except:
			pass

	return IPRdomains



#Read genes with corresponding ipr descriptions into dictionary
##########################################################################################

def IPRdescription_dictionary(IPRscan_out):
	IPRdescr = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'IPR' in entry[11]:
				if entry[0] not in IPRdescr:
					IPRdescr[entry[0]] = []
					IPRdescr[entry[0]].append(entry[12])
				else:
					if entry[12] not in IPRdescr[entry[0]]:
						IPRdescr[entry[0]].append(entry[12])
		except:
			pass

	return IPRdescr



#Read genes with corresponding pfam domains into dictionary
##########################################################################################

def PFAMdomain_dictionary(IPRscan_out):
	PFAMdomains = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'PF' in entry[4]:
				if entry[0] not in PFAMdomains:
					PFAMdomains[entry[0]] = []
					PFAMdomains[entry[0]].append(entry[4])
				else:
					if entry[4] not in PFAMdomains[entry[0]]:
						PFAMdomains[entry[0]].append(entry[4])
		except:
			pass

	return PFAMdomains



#Read genes with corresponding PFAM descriptions into dictionary
##########################################################################################

def PFAMdescription_dictionary(IPRscan_out):
	PFAMdescr = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'PF' in entry[4]:
				if entry[0] not in PFAMdescr:
					PFAMdescr[entry[0]] = []
					PFAMdescr[entry[0]].append(entry[5])
				else:
					if entry[5] not in PFAMdescr[entry[0]]:
						PFAMdescr[entry[0]].append(entry[5])
		except:
			pass

	return PFAMdescr


#Read genes with corresponding GOs (as list) into dictionary
##########################################################################################

def GO_dict(IPRscan_out):
	GOdict = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'GO' in entry[13]:
				if not entry[0] in GOdict:
					GOdict[entry[0]] = []
					GOdict[entry[0]].append(entry[13].strip())
				else:
					if entry[13] not in GOdict[entry[0]]:
						GOdict[entry[0]].append(entry[13].strip())
		except:
			pass

	return GOdict






#Read genes with corresponding EC numbers (as list) into dictionary
##########################################################################################

def ECnumber_dict(IPRscan_out):
	ECdict = {}

	for line in IPRscan_out:
		try:
			entry = line.split('\t')
			if 'KEGG' in entry[14]:
				if not entry[0] in ECdict:
					ECdict[entry[0]] = []
					ECdict[entry[0]].append(entry[14].strip())
				else:
					if entry[14] not in ECdict[entry[0]]:
						ECdict[entry[0]].append(entry[14].strip())
		except:
			pass

	return ECdict



##########################################################################################
#Parse all information into final table:
#GeneID	IPRdomain	IPRdescription	PFAMdomain	PFAMdescription	GOresults KEGGresults
##########################################################################################

#GeneList
geneIDs = make_list_geneID(infile2)

#Dictionaries
IPRdomains = IPRdomain_dictionary(infile)
IPRdescr = IPRdescription_dictionary(infile)
PFAMdomains = PFAMdomain_dictionary(infile)
PFAMdescr = PFAMdescription_dictionary(infile)
GOdict = GO_dict(infile)
ECdict = ECnumber_dict(infile)


#print(geneIDs)
#print(IPRdomains)
#print(IPRdescr)
#print(PFAMdomains)
#print(PFAMdescr)
#print(ECdict)


for geneID in geneIDs:
	if geneID in IPRdomains:
		IPRdomain = '; '.join(IPRdomains[geneID])
	else:
		IPRdomain = "-"

	if geneID in PFAMdomains:
		PFAMdomain = '; '.join(PFAMdomains[geneID])
	else:
		PFAMdomain = "-"

	if geneID in GOdict:
		GOnumber = '; '.join(GOdict[geneID])
	else:
		GOnumber = "-"

	if geneID in ECdict:
		ECnumber = '; '.join(ECdict[geneID])
	else:
		ECnumber = "-"

	if geneID in IPRdescr:
		IPRdes = '; '.join(IPRdescr[geneID])
	else:
		IPRdes = "-"

	if geneID in PFAMdescr:
		PFAMdes = '; '.join(PFAMdescr[geneID])
	else:
		PFAMdes = "-"

	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID, IPRdomain, IPRdes, PFAMdomain, PFAMdes, GOnumber, ECnumber))



outfile.close()

