#!/usr/bin/env python
#Author: Nina Dombrowski

import sys
from Bio import SeqIO

def calculate_gc(sequence):
    """Calculate the GC content of a sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return float(gc_count) / len(sequence)

def print_usage():
    """Print usage information for the script."""
    print('Usage: python3 length_gc.py example.fasta')

if len(sys.argv) != 2:
    print_usage()
    sys.exit(1)

input_file = sys.argv[1]

try:
    with open(input_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence = str(record.seq)
            gc_content = calculate_gc(sequence)
            sequence_length = len(sequence)
            #py3 code
            #print(f'{record.id}\t{gc_content:.3f}\t{sequence_length}')
            print('{}\t{:.3f}\t{}'.format(record.id, gc_content, sequence_length))
except IOError:
    #py3 code
    #print(f'Error: could not find file {input_file}')
    print('Error: could not find file {}'.format(input_file))
