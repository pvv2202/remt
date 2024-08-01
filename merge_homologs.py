'''
Merge homologs downloads from genome_from_protein.py
'''

import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Restriction enzyme and methyltransferase directories')
parser.add_argument('-o', type=str, default=None, required=True, help='Output directory name')
args = parser.parse_args()

added = []

os.mkdir(args.o)

# RE directory
dir = args.i
for dir in args.i:
    for homolog_dir in os.listdir(dir):
        if not os.path.isdir(os.path.join(dir, homolog_dir)):
            continue

        for genome in os.listdir(os.path.join(dir, homolog_dir)):
            if "annotated" not in genome or genome in added:
                continue
            shutil.copy(os.path.join(dir, homolog_dir, genome), os.path.join(args.o, genome))
            added.append(genome)
