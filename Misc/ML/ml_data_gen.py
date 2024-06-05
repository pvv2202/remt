import json
import argparse
import pandas as pd
import numpy as np
from itertools import product
import re
import os

parser = argparse.ArgumentParser(description='Generate HTML table from TSV data.')
parser.add_argument('-i', type=str, default=None, required=True, help='File input (TSV)')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
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
        self.range = 0

class Methyl(Enzyme):
    pass

class RE(Enzyme):
    pass

def calculate_similarity(string1, string2):
    set1 = set(string1)
    set2 = set(string2)

    intersection = set1.intersection(set2)
    union = set1.union(set2)

    if len(union) == 0:
        return 0

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
    if item.get("score") >= 20:
        #Filtering enzymes I don't want
        if all(keyword not in desc for keyword in ["control protein", "homing endonuclease", "subunit", "helicase", "nicking endonuclease", "orphan", "methyl-directed", "enzyme/methyltransferase"]) and "RecSeq:" in desc:
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
count = 0

for c in contigs:
    methyls = contigs[c].mts
    r_enzs = contigs[c].res
    #For each restriction enzyme and methyltransferase
    for res in r_enzs:
        r = r_enzs[res]
        first = True
        for rseq in r.seq:
            for mts in methyls:
                m = methyls[mts]
                if first:
                    for x in r_enzs:
                        if distance(min(r_enzs[x].start,r_enzs[x].end), max(r_enzs[x].start,r_enzs[x].end), min(m.start,m.end), max(m.start,m.end)) < 5000: 
                            m.range += 1
                            r_enzs[x].range += 1
                    first = False
                for mseq in m.seq:
                    target = 0
                    #If mseq <= rseq (if greater cannot block every site), run smith-waterman
                    if len(mseq) <= len(rseq):
                        aligned_mseq, aligned_rseq, score = smith_waterman(mseq, rseq)
                        if len(aligned_mseq) > 0:
                            score = score//(len(aligned_mseq))
                        else:
                            score = 0
                        asim = calculate_similarity(aligned_mseq, aligned_rseq)
                        weighted_score = round(calculate_weighted_average([m.score, r.score]) * asim)
                        if m.score >= 50 and r.score >= 50 and len(aligned_mseq) > 3 and len(aligned_rseq) > 3:
                            if (weighted_score <= 20 and mseq == rseq) or (weighted_score >= 20 and weighted_score <= 40) or weighted_score > 40:
                                target = 1
                        
                        temp = {
                            "m_score": m.score,
                            "r_score": r.score,
                            "m_seq": aligned_mseq,
                            "r_seq": aligned_rseq,
                            "smith_waterman": score,
                            "enzymes_nearby": m.range + r.range,
                            "target": target
                        }

                        hits[count] = temp

                        count += 1

df = pd.DataFrame.from_dict(hits, orient='index')

df.to_csv(args.o)