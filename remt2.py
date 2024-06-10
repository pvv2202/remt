import json
import argparse
import domainator
import pandas as pd
import numpy as np
from itertools import product
from domainator.utils import parse_seqfiles, DomainatorCDS
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import os

#TODO: PLOT GRAPH OF DISTANCE OF HITS

parser = argparse.ArgumentParser(description='Generate HTML table from JSON or TSV data.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Annotated genbank file input')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
parser.add_argument('--coords', default=True, required=False, action='store_true', help='Include coords')
parser.add_argument('--translation', default=True, required=False, action='store_true', help='Include protein sequence')
parser.add_argument('--paper', default=True, required=False, action='store_true', help='Include title of paper sequence is from')
parser.add_argument('--stats', default=None, required=False, action='store_true', help='Print stats to terminal')
parser.add_argument('--histogram', nargs='?', default=None, const=10, type=int, help='Output histogram of system distances. Specify bin size after. Default is 10')
parser.add_argument('--kde', default=None, required=False, action='store_true', help='Output kde of system distances (clipped to display positive values only)')
parser.add_argument('--excel', default=None, required=False, action='store_true', help='Make results excel friendly (default is HTML friendly)')
parser.add_argument('--filter_n', default=None, required=False, action='store_true', help='Remove hits with a lot of Ns')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
parser.add_argument('--min_distance', default=None, required=False, type=int, help='Filter by distance')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--in_range', default=5000, required=False, type=int, help='Include a metric of how many mts/res are in range of each hit')
args = parser.parse_args()

class Contig:
    def __init__(self, ref, length, topology, paper):
        self.ref = ref
        self.length = length
        self.topology = topology
        self.paper = paper
        self.fusions = 0
        self.mts = {}
        self.res = {}
        self.hits = {}

class Enzyme:
    def __init__(self, seq, start, end, score, domain, enzyme, ref, locus, translation):
        # Reference number for faster lookup
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.domain = domain
        self.enzyme = enzyme
        self.ref = ref
        self.locus = locus
        self.translation = translation
        self.range = list()
        self.matches = list()

class Methyl(Enzyme):
    pass

class RE(Enzyme):
    pass

#Storing prefixes, iupac, primes
prefixes = {"M.", "M1.", "M2.", "M3.", "M4."}

iupac = {
    'A': '[A]',
    'C': '[C]',
    'G': '[G]',
    'T': '[T]',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGTSYK]',
    'D': '[AGTWRK]',
    'H': '[ACTMWY]',
    'V': '[ACGMRS]',
    'N': '[ACGTVHDBMKWSR]'
}

def calculate_similarity(string1, string2):
    set1 = set(string1)
    set2 = set(string2)

    intersection = set1.intersection(set2)
    union = set1.union(set2)

    return len(intersection) / len(union)

def distance(A_start, A_end, B_start, B_end, topology, contig_length):
    circular_distance = 0

    if topology == "circular":
        circular_distance = contig_length - (max(A_end, B_end) - min(A_start, B_start)) # Total length - distance between start and end

    return min(max(0, B_start - A_end, A_start - B_end), circular_distance) # Return max b/c one of the start/end differences will always be negative. Then if circ dist is less return that

def iupac_match(char1, char2):
    if re.search(iupac[char1], char2):
        return True
    return False

def smith_waterman(sequence1, sequence2, exact_match_score = 3, match_score = 1, mismatch_score = -1, gap_penalty = -2):
    #Initialize the score matrix and traceback matrix
    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    score_matrix = np.zeros((rows, cols))
    traceback_matrix = np.zeros((rows, cols), dtype=int)

    #Fill the score and traceback matrices
    for i, j in product(range(1, rows), range(1, cols)):
        if sequence1[i - 1] == sequence2[j - 1]:
            to_add = exact_match_score
        elif iupac_match(sequence1[i - 1], sequence2[j - 1]):
            to_add = match_score
        else:
            to_add = mismatch_score
        match = score_matrix[i - 1, j - 1] + to_add
        delete = score_matrix[i - 1, j] + gap_penalty
        insert = score_matrix[i, j - 1] + gap_penalty
        score_matrix[i, j] = max(0, match, delete, insert)

        if score_matrix[i, j] == match:
            traceback_matrix[i, j] = 1  # Diagonal
        elif score_matrix[i, j] == delete:
            traceback_matrix[i, j] = 2  # Up
        elif score_matrix[i, j] == insert:
            traceback_matrix[i, j] = 3  # Left

    #Find the cell with the maximum score in the score matrix
    max_score = np.max(score_matrix)
    max_i, max_j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)

    #Traceback to find the alignment
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_i, max_j

    while i > 0 and j > 0 and score_matrix[i, j] > 0:
        if traceback_matrix[i, j] == 1:  #Diagonal
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 2:  #Up, gap in sequence 2
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif traceback_matrix[i, j] == 3:  #Left, gap in sequence 1
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, max_score

contigs = {}

# Populating contigs from input file
if args.i is not None:
    print("Reading Input File")
    records = parse_seqfiles(args.i)
    print("Populating Contigs")
    for i, rec in enumerate(records): # Each rec is a contig
        print(f"Reading Contig {i}", end='\r')
        count = 0
        cds_features = DomainatorCDS.list_from_contig(rec)
        for feature in cds_features:
            if len(feature.domain_features) == 0: # Skip if domainator has no annotations
                continue

            domain_features = feature.domain_features[0]

            contig_id = rec.id
            if contig_id not in contigs:
                contigs[contig_id] = (
                    Contig(
                        ref=contig_id,
                        length=len(rec.seq),
                        topology=rec.annotations['topology'],
                        paper=rec.annotations['references'][0].title
                    )
                )

            desc = domain_features.qualifiers['description'][0]
            score = domain_features.qualifiers['score'][0]

            if args.min_score is None or score >= args.min_score: # Filter by score
                # Filtering out to only consider applicable enzymes
                if args.stats and "enzyme/methyltransferase" in desc:
                    contigs[contig_id].fusions += 1
                elif all(keyword not in desc for keyword in ["control protein", "homing endonuclease", "subunit", "helicase", "nicking endonuclease", "orphan", "methyl-directed", "enzyme/methyltransferase"]) and "RecSeq:" in desc:
                    # Extracting the RecSeq
                    recseq_start = desc.find("RecSeq:") + len("RecSeq:")
                    recseq_end = desc.find(";", recseq_start)
                    seq = desc[recseq_start:recseq_end].split(", ")
                    for s in seq:
                        if "-" in s:
                            seq.remove(s)
                            continue
                        if args.filter_n:
                            if "NNN" in s:
                                seq.remove(s)
                                continue
                        if args.ignore_nonspec:
                            if len(s) < 3:
                                seq.remove(s)

                    # Extracting the EnzType
                    if "EnzType:" in desc and len(seq) > 0:
                        enz_start = desc.find("EnzType:") + len("EnzType:")
                        enz_end = desc.find(";", enz_start)
                        enz_type = desc[enz_start:enz_end]

                        # If MT/RE, respond accordingly - count is used as the reference number
                        name = domain_features.qualifiers['name'][0]
                        if any(name.startswith(prefix) for prefix in prefixes):
                            contigs[contig_id].mts[count] = (
                                Methyl(
                                   seq=seq,
                                   start=domain_features.location.start,
                                   end=domain_features.location.end,
                                   score=score,
                                   domain=name,
                                   enzyme=enz_type,
                                   ref=count,
                                   locus=domain_features.qualifiers['cds_id'][0],
                                   translation=feature.feature.qualifiers['translation'][0]
                                   )
                            )
                        else:
                            contigs[contig_id].res[count] = (
                                RE(
                                   seq=seq,
                                   start=domain_features.location.start,
                                   end=domain_features.location.end,
                                   score=score,
                                   domain=name,
                                   enzyme=enz_type,
                                   ref=count,
                                   locus=domain_features.qualifiers['cds_id'][0],
                                   translation=feature.feature.qualifiers['translation'][0]
                                   )
                            )
                        count += 1
else:
    print('Error: No input file specified.')
    exit(1)

print("\r")
print("Finding Systems")
for i, c in enumerate(contigs.values()):
    print(f'{int(i / len(contigs) * 100)}%', end='\r')
    # For each restriction enzyme and methyltransferase
    for r in c.res.values():
        for rseq in r.seq:
            for m in c.mts.values():
                for mseq in m.seq:
                    dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length)
                    # Add any re/mts within specified range
                    if args.in_range and dist < args.in_range:
                        m.range.append(f'{r.domain}: {r.seq}')
                        r.range.append(f'{m.domain}: {m.seq}')

                    # If mseq <= rseq (if greater cannot block every site), run smith-waterman
                    if len(mseq) <= len(rseq):
                        aligned_mseq, aligned_rseq, score = smith_waterman(mseq, rseq)
                        #Filtering by score, length of rseq, absolute value of mseq
                        if score > 2 and len(aligned_rseq) > 3 and (len(aligned_mseq) == len(mseq)):
                            sim = calculate_similarity(mseq, rseq)
                            asim = calculate_similarity(aligned_mseq, aligned_rseq)
                            if (mseq != rseq and "NNN" not in mseq and "NNN" not in rseq) or mseq == rseq:
                                #Doing a better job getting the score
                                weighted_score =  round((asim + sim) * score) + round(1000/dist) if dist > 0 else 1000

                                if asim == 1 and r.score > 80 and m.score > 80:
                                    if r.ref not in c.hits or weighted_score > c.hits[r.ref]["Score"]:
                                        temp = {
                                            "Contig": c.ref,
                                            "Types": [m.enzyme, r.enzyme],
                                            "Domains": [m.domain, r.domain],
                                            "Sequences": [mseq, rseq],
                                            "Alignment": [aligned_mseq, aligned_rseq],
                                            "Loci": [m.locus, r.locus],
                                            "Score": weighted_score,
                                            "Scores": [m.score, r.score],
                                            "Distance": dist,
                                        }

                                        if args.in_range:
                                            if len(r.range) == 0:
                                                r.range.append("None")
                                            if len(m.range) == 0:
                                                m.range.append("None")
                                            for mt_in_range in r.range:
                                                temp["MT within " + str(args.in_range) + " bp"] = mt_in_range
                                            for rt_in_range in m.range:
                                                temp["RE within " + str(args.in_range) + " bp"] = rt_in_range

                                        if args.coords:
                                            temp["Coords"] = [str(m.start) + ", " + str(m.end), str(r.start) + ", " + str(r.end)]

                                        if args.translation:
                                            temp["Translation"] = [m.translation, r.translation]

                                        if args.paper:
                                            temp["Paper"] = c.paper

                                        c.hits[r.ref] = temp

print("\r    ")

hits = {}
for c in contigs:
    hits.update(contigs[c].hits)

#Filtering hits for distance
if args.min_distance is not None:
    hits = {h: v for h, v in hits.items() if v["Distance"] >= args.min_distance}

#Printing necessary data for stats
if args.stats:
    num_none = 0
    num_both = 0
    num_only_re = 0
    num_only_mt = 0
    hit_stats = []
    re_only = []

    for c in contigs.values():
        if (len(c.res) > 0 and len(c.mts) > 0) or c.fusions > 0:
            num_both += 1
            curr = c.fusions + len(c.hits)
            if len(hit_stats) < curr + 1:
                hit_stats.extend([0] * (curr - len(hit_stats) + 1))
            hit_stats[curr] += 1
        elif len(c.res) == 0 and len(c.mts) == 0:
            num_none += 1
        elif len(c.res) > 0 and len(c.mts) == 0:
            num_only_re += 1
            if c not in re_only:
                re_only.append(c)
        elif len(c.res) == 0 and len(c.mts) > 0:
            num_only_mt += 1

    print("Data for Contigs with both MT and RE:")
    for i, h in enumerate(hit_stats):
        print(f"Contigs with {i} Systems: {h}")

    print("\nData for Contigs in general:")
    print("Num Contigs: " + str(len(contigs)))
    print("Num Both/Fusion: " + str(num_both))
    print("Num None: " + str(num_none))
    print("Num Only RE: " + str(num_only_re))
    print("Num Only MT: " + str(num_only_mt))

if args.histogram is not None or args.kde:
    distances = [hit["Distance"] for hit in hits.values()]
    if args.histogram:
        # Plot histogram
        plt.hist(distances, bins=args.histogram, alpha=0.75, edgecolor='black', density=False)
        plt.title('System Distances')
        plt.xlabel('Distance')
        plt.ylabel('Count')
        plt.show()
    if args.kde:
        # Plot KDE
        sns.kdeplot(distances, fill=True, clip=(0, None))
        plt.title('KDE of Distances')
        plt.xlabel('Distance')
        plt.ylabel('Density')
        plt.show()

if hits:
    # Get the keys from the first dictionary in hits
    for key, value in hits.items():
        if isinstance(value, dict):
            header_keys = value.keys()
            break
    else:
        # Handle the case when no dictionaries are found in hits
        header_keys = []

    # Generate the HTML table with centered content
    html_table = "<table id='myTable' style='border-collapse: collapse; font-family: Courier New, monospace; font-size: 0.7em'>\n"
    # Create table header
    html_table += "<thead><tr>"
    for key in header_keys:
        html_table += "<th style='border: 1px solid black; padding: 8px; text-align: left;'>" + str(key) + "</th>"
    html_table += "</tr></thead>\n"
    # Create table rows
    html_table += "<tbody>"
    for item in hits.values():
        if isinstance(item, dict):
            html_table += "<tr>"
            for key, value in item.items():
                if key == "Alignment" or key == "Sequences" or key == "Coords" or key == "Scores" or key == "Types" or key == "Domains" or key == "Loci" or key == "Translation":
                    if args.excel:
                        html_table += "<td style='border: 1px solid black; padding: 8px; text-align: left;'>" + "R: " + str(value[1]) + "&#10;" + "M: " + str(value[0]) + "</td>"
                    else:
                        html_table += "<td style='border: 1px solid black; padding: 8px; text-align: left;'>" + "R: " + str(value[1]) + "<br>" + "M: " + str(value[0]) + "</td>"
                else:
                    html_table += "<td style='border: 1px solid black; padding: 8px; text-align: left;'>" + str(value) + "</td>"
            html_table += "</tr>\n"
    html_table += "</tbody></table>"

    output_filename = args.o + ".html"

    # Check if the file already exists
    if os.path.exists(output_filename):
        mode = "w"  # If the file exists, open it in write mode to overwrite
    else:
        mode = "a"  # If the file doesn't exist, open it in append mode

    with open(output_filename, mode) as file:
        file.write(html_table)

    html_script = """
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.3/css/jquery.dataTables.min.css">
    <script src="https://cdn.datatables.net/1.11.3/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready(function() {
            $('#myTable').DataTable();
        });
    </script>
    """
    # Append the JavaScript code to the HTML file
    with open(output_filename, "a") as file:
        file.write(html_script)