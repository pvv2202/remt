'''
Sort through homologs of an MT and RE to output distances of systems and whether other enzymes are nearby
'''

import argparse
import os
from utils import Annotations, distance

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Restriction enzyme and methyltransferase directories')
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
            for m in c.mts.values():
                if distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length) < 10000:
                    fasta += f">c__{contig_id}_{r.enzyme}_{genome}\n{r.translation}\n"
                else:
                    fasta += f">s__{contig_id}_{r.enzyme}_{genome}\n{r.translation}\n"

# Writing to genomes file
with open(f'{args.o}_genomes', 'w') as file:
    file.write(genomes)

# Writing to fasta files
with open(f'{args.o}.fasta', 'w') as file:
    file.write(fasta)
