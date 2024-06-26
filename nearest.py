import json
import argparse
import pandas as pd
import numpy as np
from numba import jit
from itertools import product
import re
import os

parser = argparse.ArgumentParser(description='Generate HTML table from JSON or TSV data.')
parser.add_argument('-i', type=str, default=None, required=True, help='File input (JSON or TSV)')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (JSON and HTML)')
parser.add_argument('--coords', default=None, required=False, action='store_true', help='Include coords')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
args = parser.parse_args()

# Handling input file
if args.i is not None:
    in_file = args.i
    if in_file.endswith('.tsv'):
        data = pd.read_csv(in_file, sep='\t')
    else:
        print('Error: Unsupported file format. Please provide a TSV file.')
        exit(1)
else:
    print('Error: No input file specified.')
    exit(1)

#Storing prefixes, iupac, primes
prefixes = {"M.", "M1.", "M2.", "M3.", "M4.", "M5.", "M6."}

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

class Contig:
    def __init__(self, ref):
        self.ref = ref
        self.mts = []
        self.res = []

class Enzyme:
    def __init__(self, seq, start, end, score, domain, locus, ref):
        # Reference number for faster lookup
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.domain = domain
        self.locus = locus
        self.ref = ref
        self.range = list()

class Methyl(Enzyme):
    pass

class RE(Enzyme):
    pass

def calculate_similarity(string1, string2):
    set1 = set(string1)
    set2 = set(string2)

    intersection = set1.intersection(set2)
    union = set1.union(set2)

    return len(intersection) / len(union)

def calculate_weighted_average(scores):
    weights = [(1 / (score - min(scores) + 1)) for score in scores]
    weighted_sum = sum(score * weight for score, weight in zip(scores, weights))
    total_weight = sum(weights)

    return weighted_sum / total_weight

def distance(A_start, A_end, B_start, B_end):
    #Circularity = min(A_start, B_start) + contig_length - max(A_end, B_end)
    return max(0, B_start - A_end, A_start - B_end)

def iupac_match(char1, char2):
    if re.search(iupac[char1], char2):
        return True
    return False

def smith_waterman(sequence1, sequence2, exact_match_score = 3, match_score = 1, mismatch_score = -1, gap_penalty = -2):
    # Initialize the score matrix and traceback matrix
    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    score_matrix = np.zeros((rows, cols))
    traceback_matrix = np.zeros((rows, cols), dtype=int)

    # Fill the score and traceback matrices
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

    # Find the cell with the maximum score in the score matrix
    max_score = np.max(score_matrix)
    max_i, max_j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)

    # Traceback to find the alignment
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_i, max_j

    while i > 0 and j > 0 and score_matrix[i, j] > 0:
        if traceback_matrix[i, j] == 1:  # Diagonal
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 2:  # Up, gap in sequence 2
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif traceback_matrix[i, j] == 3:  # Left, gap in sequence 1
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, max_score

contigs = {}

count = 0
#Populate contigs by iterating over data. Now using iterrows() to speed up (count is the reference used in dictionary)
for index, item in data.iterrows():
    desc = item.get("domain_descriptions")
    ref = item.get("contig")
    if args.min_score is None or item.get("score") >= args.min_score:
        #Filtering enzymes I don't want
            if ref not in contigs:
                contigs[ref] = Contig(ref)
            seq = desc
            #If MT/RE, respond accordingly - count is used as the reference number
            if any(item.get("domain").startswith(prefix) for prefix in prefixes):
                contigs[ref].mts.append(Methyl(seq, min(item.get("start"), item.get("end")), max(item.get("start"), item.get("end")), item.get("score"), item.get("domain"), item.get("cds"), count))
                contigs[ref].mts.sort(key = lambda x: x.start)
            else:
                contigs[ref].res.append(RE(seq, min(item.get("start"), item.get("end")), max(item.get("start"), item.get("end")), item.get("score"), item.get("domain"), item.get("cds"), count))
                contigs[ref].res.sort(key = lambda x: x.start)
            count += 1

hits = {}

for c in contigs:
    m = contigs[c].mts
    r = contigs[c].res
    mPrev = None
    j = 0
    for i, rz in enumerate(r):
        while j < len(m) - 1 and m[j+1].start < rz.start:
            j += 1
        if j < len(m) - 1 and rz.start - m[j].end < m[j+1].start - rz.end:
            rToAdd = rz
            if mPrev is not None and m[j].start == mPrev.start:
                if rz.start - mPrev.end > mPrev.start - r[i-1].start:
                    rToAdd = r[i-1]
            mToAdd = m[j]
            dist = distance(rToAdd.start, rToAdd.end, mToAdd.start, mToAdd.end)
            if dist > 5000 and calculate_similarity(rToAdd.domain, mToAdd.domain) > 0.8:
                temp = {
                    "Contig": c,
                    "Loci": [mToAdd.locus, rToAdd.locus],
                    "Domains": [mToAdd.domain, rToAdd.domain],
                    "Sequences": [mToAdd.seq, rToAdd.seq],
                    "Scores": [mToAdd.score, rToAdd.score],
                    "Coords": [str(mToAdd.start) + ", " + str(mToAdd.end), str(rToAdd.start) + ", " + str(rToAdd.end)],
                    "Distance": dist
                }
                hits[rToAdd.ref] = temp
        else:
            mPrev = m[j+1] if j < len(m) - 1 else None



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
                if key == "Alignment" or key == "Sequences" or key == "Coords" or key == "Scores" or key == "Types" or key == "Domains" or key == "Loci":
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