"""
Sort through homologs of an MT and RE to output distances of systems and whether other enzymes are nearby
"""

import argparse
import os
from utils import Annotations, distance

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Merged RE/MT directory')
parser.add_argument('--type', type=str, default="RE", required=True, help='Enzyme type')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
parser.add_argument('-o', type=str, default=None, required=True, help='Output file name')
args = parser.parse_args()

fasta = ""
genomes = ""

# RE directory
dir = args.i
for genome in os.listdir(dir):
    genomes += f"{genome}\n"

    # Reformat args to be passed to utils, then parse
    temp_args_dict = vars(args).copy()
    temp_args_dict['i'] = os.path.join(dir, genome)
    temp_args = argparse.Namespace(**temp_args_dict)

    annotation = Annotations(temp_args)
    contigs = annotation.parse()

    for contig_id, c in contigs.items():
        for r in c.res.values():
            split = None
            found_close = False
            for m in c.mts.values():
                if distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length) < 10000:
                    if args.type == "RE":
                        fasta += f">c__{contig_id}_{r.enzyme}_{genome}\n{r.translation}\n"
                        found_close = True
                        break
                    elif args.type == "MT":
                        fasta += f">c__{contig_id}_{m.enzyme}_{genome}\n{m.translation}\n"
                        found_close = True
                        break
                    else:
                        raise Exception(f"Unknown type {args.type}. Expects RE or MT")
                else:
                    if args.type == "RE":
                        split = f">s__{contig_id}_{r.enzyme}_{genome}\n{r.translation}\n"
                    elif args.type == "MT":
                        split = f">s__{contig_id}_{m.enzyme}_{genome}\n{m.translation}\n"
                    else:
                        raise Exception(f"Unknown type {args.type}. Expects RE or MT")

            if split and not found_close:
                fasta += split

# Writing to genomes file
with open(f'{args.o}_genomes', 'w') as file:
    file.write(genomes)

# Writing to fasta files
with open(f'{args.o}.fasta', 'w') as file:
    file.write(fasta)
