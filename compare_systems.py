"""
Program that finds close systems and matches them to homologous split systems
"""

import argparse
import pickle
import os
import csv
from domainator.seq_dist import seq_dist
from domainator.data_matrix import DataMatrix
from utils import Annotations, distance

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Annotated genbank or pickle file input')
parser.add_argument('--mats', nargs='+', type=str, default=None, required=False, help='Similarity matrices. RE goes first')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
args = parser.parse_args()

contigs = {}
re_sequences=""
mt_sequences=""

if args.i.endswith('.gb'):
    Annotations = Annotations(args)
    contigs = Annotations.parse()

    # Save the contigs dictionary as a pickle object
    with open(f'{args.o}.pkl', 'wb') as f:
        pickle.dump(contigs, f)
elif args.i.endswith(".pkl"):
    with open(f'{args.i}', 'rb') as f:
        contigs = pickle.load(f)
else:
    print('Error: Unsupported file format. Please provide a .gb or .pkl file.')
    exit(1)

os.makedirs(args.o, exist_ok=True)

if not args.mats: # Only do this if the matrices are not provided
    for i, contig in enumerate(contigs.values()):
        print(f'{int(i / len(contigs) * 100)}%', end='\r')
        for r in contig.res.values():
            re_sequences += f">{r.locus}\n{r.translation}\n"
        for m in contig.mts.values():
            mt_sequences += f">{m.locus}\n{m.translation}\n"

    # Write fasta files, generate matrices, load matrices
    with open(f"{args.o}/{args.o}_re_sequences.fasta", "w") as f:
        f.write(re_sequences)
    with open(f"{args.o}/{args.o}_mt_sequences.fasta", "w") as f:
        f.write(mt_sequences)

    seq_dist(f"{args.o}/{args.o}_re_sequences.fasta", "fasta", f"{args.o}/{args.o}_re_sequences.fasta", "fasta", None, "diamond_us", "score", 8, None, None, f"{args.o}/{args.o}_re_similarity_matrix.hdf5", 0)
    seq_dist(f"{args.o}/{args.o}_mt_sequences.fasta", "fasta", f"{args.o}/{args.o}_mt_sequences.fasta", "fasta", None, "diamond_us", "score", 8, None, None, f"{args.o}/{args.o}_mt_similarity_matrix.hdf5", 0)

re_mat_path = f"{args.o}/{args.o}_re_similarity_matrix.hdf5"
mt_mat_path = f"{args.o}/{args.o}_mt_similarity_matrix.hdf5"

if args.mats:
    # RE first, MT second
    re_mat_path = args.mats[0]
    mt_mat_path = args.mats[1]

re_matrix = DataMatrix.from_file(re_mat_path)
mt_matrix = DataMatrix.from_file(mt_mat_path)

all_re = {r.locus: (c.ref, r) for c in contigs.values() for r in c.res.values()}
hits = {}

# Compare the systems
count = 0
for i, c in enumerate(contigs.values()):
    print(f'{int(i / len(contigs) * 100)}%', end='\r')
    for r in c.res.values():
        for m in c.mts.values():
            # If distance under threshold, check for matching re and then for matching mt within the right distance. If all applies, add to hits
            s1_dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length)
            if s1_dist > 2500:
                continue
            for r2_loc in all_re:
                if r.locus not in re_matrix.rows or r2_loc not in re_matrix.columns:
                    continue
                if re_matrix.data[re_matrix.row_to_idx[r.locus], re_matrix.column_to_idx[r2_loc]] < 50:
                    continue
                r2 = all_re[r2_loc][1]
                c2 = contigs[all_re[r2_loc][0]]
                for m2 in c2.mts.values(): # Check mts on the re's contig
                    if m.locus not in mt_matrix.rows or m2.locus not in mt_matrix.columns:
                        continue
                    if mt_matrix.data[mt_matrix.row_to_idx[m.locus], mt_matrix.column_to_idx[m2.locus]] < 50:
                        continue
                    s2_dist = distance(min(r2.start, r2.end), max(r2.start, r2.end), min(m2.start, m2.end), max(m2.start, m2.end), c2.topology, c2.length)
                    if s2_dist < 5000:
                        break # If we find an MT nearby, this RE is not in a split system
                    hits[count] = {
                        "Contig": f"M1: {c.ref}\nR1: {c2.ref}",
                        "Strain": f"M1: {c.strain}\nR1: {c2.strain}",
                        "Types": f"M1: {m.enzyme}\nR1: {r.enzyme}\nM2: {m2.enzyme}\nR2: {r2.enzyme}",
                        "Domains": f"M1: {m.domain}\nR1: {r.domain}\nM2: {m2.domain}\nR2: {r2.domain}",
                        "Distance": f"M1: {s1_dist}\nR1: {s2_dist}",
                        "Sequences": f"M1: {m.seq}\nR1: {r.seq}\nM2: {m2.seq}\nR2: {r2.seq}",
                        "Loci": f"M1: {m.locus}\nR1: {r.locus}\nM2: {m2.locus}\nR2: {r2.locus}",
                        "Scores": f"M1: {m.score}\nR1: {r.score}\nM2: {m2.score}\nR2: {r2.score}",
                        "Translation": f"M1: {m.translation}\nR1: {r.translation}\nM2: {m2.translation}\nR2: {r2.translation}"
                    }
                    count += 1

if hits:
    # Get the keys from the first dictionary in hits
    header_keys = hits[0].keys()

    output_filename = args.o + ".csv"
    output_filepath = os.path.join(args.o, output_filename)

    # Write the CSV file
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header_keys)
        writer.writeheader()

        for hit in hits.values():
            # Flatten the nested lists for CSV output
            flattened_hit = {key: "\n".join(map(str, value)) if isinstance(value, list) else value for key, value in hit.items()}
            writer.writerow(flattened_hit)

    print(f"CSV file saved to {output_filepath}")