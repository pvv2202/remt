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
parser.add_argument('--stats', default=None, required=False, action='store_true', help='Print stats to terminal')
parser.add_argument('--filter_n', default=None, required=False, action='store_true', help='Remove hits with a lot of Ns')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
parser.add_argument('--min_distance', default=None, required=False, type=int, help='Filter by distance')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--in_range', default=None, required=False, type=int, help='Include a metric of how many mts/res are in range of each hit')
args = parser.parse_args()

# Handling input file
if args.i is not None:
    in_file = args.i
    if in_file.endswith('.json'):
        with open(in_file, 'r') as file:
            data = json.load(file)
    elif in_file.endswith('.tsv'):
        data = pd.read_csv(in_file, sep='\t')
    else:
        print('Error: Unsupported file format. Please provide a JSON or TSV file.')
        exit(1)
else:
    print('Error: No input file specified.')
    exit(1)

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

class Contig:
    def __init__(self, ref):
        self.ref = ref
        self.hits = 0
        self.fusions = 0
        self.mts = {}
        self.res = {}

class Enzyme:
    def __init__(self, seq, start, end, score, domain, enzyme, ref, locus):
        # Reference number for faster lookup
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.domain = domain
        self.enzyme = enzyme
        self.ref = ref
        self.locus = locus
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
        #Factoring enzyme/methyltransferase for the stats
        if args.stats and "enzyme/methyltransferase" in desc:
            if ref not in contigs:
                contigs[ref] = Contig(ref)
            contigs[ref].fusions += 1
            contigs[ref].hits += 1
        #Filtering enzymes I don't want
        elif all(keyword not in desc for keyword in ["control protein", "homing endonuclease", "subunit", "helicase", "nicking endonuclease", "orphan", "methyl-directed", "enzyme/methyltransferase"]) and "RecSeq:" in desc:
            if ref not in contigs:
                contigs[ref] = Contig(ref)
            i = desc.find("RecSeq:")
            recseq_start = i + len("RecSeq:")
            recseq_end = desc.find(";", recseq_start)
            long = desc[recseq_start:recseq_end]
            if ", " in long:
                seq = long.split(", ")
            else:
                seq = [long]
            
            if args.filter_n:
                for s in seq:
                    if "NNN" in s or not args.ignore_nonspec or len(s) < 3:
                        seq.remove(s)

            if "EnzType:" in desc and len(seq) > 0:
                j = desc.find("EnzType:")
                et_start = j + len("EnzType:")
                et_end = desc.find(";", et_start)
                et = desc[et_start:et_end]
                #If MT/RE, respond accordingly - count is used as the reference number
                if any(item.get("domain").startswith(prefix) for prefix in prefixes):
                    contigs[ref].mts[count] = Methyl(seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain"), et, count, item.get("cds"))
                else:
                    contigs[ref].res[count] = RE(seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain"), et, count, item.get("cds"))
                count += 1

hits = {}

for c in contigs:
    methyls = contigs[c].mts
    r_enzs = contigs[c].res
    #For each restriction enzyme and methyltransferase
    for res in r_enzs:
        r = r_enzs[res]
        best = [r.seq, None, -30]
        best_match = [None, 0, None]
        for rseq in r.seq:
            for mts in methyls:
                m = methyls[mts]
                for mseq in m.seq:
                    #If args in range is not none, check if the MT is in range. If it is, 
                    if args.in_range is not None:
                        if distance(min(r.start,r.end), max(r.start,r.end), min(m.start,m.end), max(m.start,m.end)) < args.in_range: 
                            sim = calculate_similarity(rseq, mseq)
                            if best_match[0] is None:
                                best_match = [rseq, mseq, sim]
                            elif sim > best_match[2]:
                                best_match = [rseq, mseq, sim]
                    #If mseq <= rseq (if greater cannot block every site), run smith-waterman
                    if len(mseq) <= len(rseq):
                        aligned_mseq, aligned_rseq, score = smith_waterman(mseq, rseq)
                        #Filtering by score, length of rseq, absoltue value of mseq
                        if score > 2 and len(aligned_rseq) > 3 and (len(aligned_mseq) == len(mseq)):
                            dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end))
                            sim = calculate_similarity(m.seq, r.seq)
                            asim = calculate_similarity(aligned_mseq, aligned_rseq)
                            if r.ref in hits:
                                csim = calculate_similarity(hits[r.ref]["Alignment"][0], hits[r.ref]["Alignment"][1])
                            #More filters
                            if r.ref not in hits or (asim > csim - 0.2 and dist < hits[r.ref]["Distance"]):
                                if (mseq != rseq and "NNN" not in mseq and "NNN" not in rseq) or mseq == rseq:
                                    #Doing a better job getting the score
                                    weighted_score = round(calculate_weighted_average([m.score, r.score]) * sim)
                                    if (weighted_score <= 20 and mseq == rseq) or (weighted_score >= 20 and weighted_score <= 40 and sim > 0.9) or weighted_score > 40:
                                        if weighted_score > best[2]:
                                            best = [r, m, weighted_score, rseq, mseq, aligned_rseq, aligned_mseq]
                if best_match[0] is not None:
                    m.range.append(best_match[0])
        if best_match[0] is not None:
            r.range.append(best_match[1])
        if best[1] is not None:
            r = best[0]
            m = best[1]
            temp = {
                "Contig": c,
                "Types": [m.enzyme, r.enzyme],
                "Domains": [m.domain, r.domain],
                "Sequences": [best[4], best[3]],
                "Alignment": [best[6], best[5]],
                "Loci": [m.locus, r.locus],
                "Score": best[2],
                "Scores": [m.score, r.score],
                "Distance": distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end))
            }
            if args.in_range is not None:
                for mt_in_range in r.range:
                    temp["MT within " + str(args.in_range) + " bp"] = mt_in_range
                for rt_in_range in m.range:
                    temp["RE within " + str(args.in_range) + " bp"] = rt_in_range

            if args.coords:
                temp["Coords"] = [str(m.start) + ", " + str(m.end), str(r.start) + ", " + str(r.end)]

            hits[r.ref] = temp

        
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

    for contig in contigs:
        c = contigs[contig]
        for h in hits:
            if c.ref == hits[h]["Contig"]:
                c.hits += 1

        if len(hit_stats) < c.hits + 1:
            hit_stats.extend([0] * (c.hits - len(hit_stats) + 1))
        hit_stats[c.hits] += 1

        if (len(c.res) > 0 and len(c.mts) > 0) or c.fusions > 0:
            num_both += 1
        elif len(c.res) == 0 and len(c.mts) == 0:
            num_none += 1
        elif len(c.res) > 0 and len(c.mts) == 0:
            num_only_re += 1
            if c not in re_only:
                re_only.append(c)
        elif len(c.res) == 0 and len(c.mts) > 0:
            num_only_mt += 1
        #Add re/mt fusion data

    for i, h in enumerate(hit_stats):
        print(str(i) + ": " + str(h))

    print("Num Contigs: " + str(len(contigs)))
    print("Num Both/Fusion: " + str(num_both))
    print("Num None: " + str(num_none))
    print("Num Only RE: " + str(num_only_re))
    print("Num Only MT: " + str(num_only_mt))

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