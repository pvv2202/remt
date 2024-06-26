import json
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, default=None, required=True, help='File input (JSON)')
args = parser.parse_args()

if args.i is not None:
    in_file = args.i
    with open(in_file, 'r') as file:
        data = json.load(file)
else:
    print('Error: No input file specified.')
    exit(1)

gb_paths_prim = []
gb_paths_sec = []

for item in data:
    if 'gbk_fns' in item:
        gb_paths_prim.extend(item['gbk_fns'])

for item in data:
    if 'gbk_fns_sec' in item:
        gb_paths_sec.extend(item['gbk_fns_sec'])

def concatenate_genbank_files(input_files, output_file):
    command = f"cat {' '.join(input_files)} > {output_file}"
    os.system(command)

concatenate_genbank_files(gb_paths_prim, "concatenated_primary.gb")
concatenate_genbank_files(gb_paths_sec, "concatenated_secondary.gb")