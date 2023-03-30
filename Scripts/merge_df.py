#python3 script
#author: Nina Dombrowski

import pandas as pd
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(
    description='Merge multiple TSV files into a single file',
    epilog='Example: python merge_df.py */[A-Z]_*tsv accession Annotations.txt',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument('files', nargs='+', help='List of files to merge')
parser.add_argument('column', help='Name of the column to merge on')
parser.add_argument('output', help='Name of the output file')
args = parser.parse_args()

#sort the input files
file_names = sorted(args.files, key=lambda x: x.split('/')[1].split('_')[0])

# Load the first file into a pandas DataFrame
df = pd.read_csv(file_names[0], delimiter='\t')

# Iterate over the remaining files and merge their data into the DataFrame
for file in file_names[1:]:
    temp_df = pd.read_csv(file, delimiter='\t')
    df = pd.merge(df, temp_df, on=args.column, how='outer')

# Write the merged DataFrame to a new file
df.to_csv(args.output, sep='\t', index=False)
