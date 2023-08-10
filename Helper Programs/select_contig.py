import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('-i', type=str, default=None, required=True, help='File inputs (TSV), by domain first')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
args = parser.parse_args()

# Handling input file
if args.i is not None:
    if args.i.endswith('.tsv'):
        file = pd.read_csv(args.i, sep='\t')
    else:
        print('Error: Unsupported file format. Please provide a TSV file.')
        exit(1)

contigs = {}
for item in file:
    contig = item.get("contig")
    if contig not in contigs:
        contigs[contig] = True

for key in contigs:
    print(key)