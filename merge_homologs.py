"""
Merge homologs downloaded from genome_from_protein.py. Useful to get all homologs in one place.
"""

import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Restriction enzyme and methyltransferase directories')
parser.add_argument('-o', type=str, default=None, required=True, help='Output directory name')
args = parser.parse_args()

# added = []
#
# os.makedirs(args.o, exist_ok=True)
#
# def create_dir_set(dir):
#     dir_set = set()
#     for homolog_dir in os.listdir(dir): # Iterate over homologs
#         if not os.path.isdir(os.path.join(dir, homolog_dir)):
#             continue
#
#         for genome in os.listdir(os.path.join(dir, homolog_dir)): # Iterate over genomes
#             if "annotated" in genome:
#                 dir_set.add(genome)
#
#     return dir_set
#
# # Create sets
# re_set = create_dir_set(args.i[0])
# mt_set = create_dir_set(args.i[1])
#
# # Find intersect
# intersect = re_set.intersection(mt_set)
#
# # Find all the paths from RE and move them to a new folder
# for homolog_dir in os.listdir(args.i[0]):
#     if not os.path.isdir(os.path.join(args.i[0], homolog_dir)):
#         continue
#
#     for genome in os.listdir(os.path.join(args.i[0], homolog_dir)):  # Iterate over genomes
#         if genome in intersect and "annotated" in genome:
#             shutil.copy(os.path.join(args.i[0], homolog_dir, genome), os.path.join(args.o, genome))

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
