#!/usr/bin/python
"""
Author: Anja Spang
Recodes amino acids to Dayoff 4, Dayoff 6 or SR4 categories
September 2018

#Note: selenocysteine will be recoded to a gap

Dayhoff6
A,G,P,S,T = A;  D,E,N,Q = D; H,K,R = H; F,Y,W = F; I,L,M,V = I; C = C

Dayhoff4
A,G,P,S,T = A; D,E,N,Q = D; H,K,R = H; F,Y,W,I,L,M,V = F; C = ?

Sr4
A,G,N,P,S,T = A; C,H,W,Y = C; D,E,K,Q,R = G; F,I,L,M,V = T
"""

"""
D6
>test
AGPSTDENQHKRFYWILMVC
>test
AAAAADDDDHHHFFFIIIIC

S4:
>test
AGPSTDENQHKRFYWILMVC
>test
AAAAAGGAGCGGTCCTTTTC


"""

#Import modules
##########################################################################################

import sys
from Bio import SeqIO
import argparse



#ARGPARSE
##########################################################################################
parser = argparse.ArgumentParser(description='Recode amino acids to Dayoff 4, Dayoff 6 or SR4 categories')
parser.add_argument('-a','--input1', metavar='fasta', required = True, help='multifasta')
parser.add_argument('-f','--input2', metavar= str, required = True, help='category, takes three arguments; either D4, D6 or S4')
parser.add_argument('-o','--output', metavar='fasta', required = True, help= 'multifasta with replaced character')
args = parser.parse_args()

if args.input1 is None and not os.path.exists(args.input):
  print("multifasta file is file missing!")
  sys.exit(-1)

if args.input2 is None and not os.path.exists(args.input):
  print("character to be found is missing")
  sys.exit(-1)




#Global variables
##########################################################################################

infile = open(args.input1)
outfile= open("%s" % args.output,'w')
category= str(args.input2)




#Function for recoding 
##########################################################################################

def recode_seq(multifasta, cat):
    DH6 = {'A':'A','G':'A','P':'A','S':'A','T':'A','D':'D','E':'D','N':'D','Q':'D','H':'H','K':'H','R':'H','F':'F','Y':'F','W':'F','I':'I','L':'I','M':'I','V':'I','C':'C'}
    DH4 = {'A':'A','G':'A','P':'A','S':'A','T':'A','D':'D','E':'D','N':'D','Q':'D','H':'H','K':'H','R':'H','F':'F','Y':'F','W':'F','I':'F','L':'F','M':'F','V':'F','C':'-'}
    SR4 = {'A':'A','G':'A','N':'A','P':'A','S':'A','T':'A','C':'C','H':'C','W':'C','Y':'C','D':'G','E':'G','K':'G','Q':'G','R':'G','F':'T','I':'T','L':'T','M':'T','V':'T'}
    
    for seq_record in SeqIO.parse(multifasta,"fasta"):
        header = str(seq_record.description).strip()
        seq = str(seq_record.seq).strip()
        new_seq = ''
        
        for item in seq:
            if cat == "D6":
                if item in DH6:
                    new_seq = new_seq+str(DH6.get(item))
                else: 
                    new_seq = new_seq+str('-')
            elif cat == "D4":
                if item in DH4:
                    new_seq = new_seq+str(DH4.get(item))	
                else:
                    new_seq = new_seq+str('-')
            elif cat == "S4":
                if item in SR4:
                    new_seq = new_seq+str(SR4.get(item))
                else:
                    new_seq = new_seq+str('-')
				
				
        outfile.write(">%s\n%s\n" % (header, new_seq))


#Call main
##########################################################################################

def main():
	recode_seq(infile, category)


if __name__ == "__main__":
	main()




#Close files
##########################################################################################
outfile.close()
