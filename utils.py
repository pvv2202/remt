import gb_io
import numpy as np
from numba import jit

class Annotations:
    class Contig:
        def __init__(self, ref, length, topology, strain):
            self.ref = ref
            self.length = length
            self.topology = topology
            self.strain = strain
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

    def __init__(self, args):
        self.args = args
        self.prefixes = {"M.", "M1.", "M2.", "M3.", "M4."}
        self.contigs = {}

    def extract_domainator_cds(self, record):
        cdss = {}
        for feature in record.features:
            qualifiers = {q.key: q.value for q in feature.qualifiers}
            if "cds_id" not in qualifiers:
                continue
            if qualifiers["cds_id"] not in cdss:
                cdss[qualifiers["cds_id"]] = {}
            dict = cdss[qualifiers["cds_id"]]
            if feature.kind == "CDS":
                dict["translation"] = qualifiers["translation"]
            elif feature.kind == "Domainator":
                dict["start"] = feature.location.start
                dict["end"] = feature.location.end
                dict["desc"] = qualifiers["description"]
                dict["name"] = qualifiers["name"]
                dict["score"] = qualifiers["score"]

        temp = cdss.copy()
        for id, cds in cdss.items():
            if len(cds.keys()) < 2:
                del temp[id]
        return temp

    def process_record(self, record, count):
        cds_features = self.extract_domainator_cds(record)
        for cds_id, feature in cds_features.items():
            contig_id = record.name
            desc = feature["desc"].split(";")
            if contig_id not in self.contigs:
                strain = desc[3][9:]
                self.contigs[contig_id] = Annotations.Contig(
                    ref=contig_id,
                    length=record.length,
                    topology="circular" if record.circular else "linear",
                    strain=strain
                )

            if self.args.min_score and feature["score"] < self.args.min_score:
                continue

            enz_type = desc[0][7:]

            if "enzyme/methyltransferase" in enz_type:
                self.contigs[contig_id].fusions += 1
            elif all(keyword not in enz_type for keyword in ["control protein", "homing endonuclease", "subunit", "helicase", "nicking endonuclease", "orphan", "methyl-directed"]) and "RecSeq:" in desc[1]:
                seq = desc[1][7:].split(", ")
                # Filter recognition sequence. Remove those with gaps, and if specified those with a long chain of N's, or less than 3 bases
                seq = [s for s in seq if "-" not in s and (not self.args.ignore_nonspec or (len(s) >= 3 and "NNN" not in s))]

                if not seq:
                    continue

                enzyme_class = Annotations.Methyl if any(feature["name"].startswith(prefix) for prefix in self.prefixes) else Annotations.RE
                enzyme_dict = self.contigs[contig_id].mts if enzyme_class == Annotations.Methyl else self.contigs[contig_id].res
                enzyme_dict[count] = enzyme_class(
                    seq=seq,
                    start=feature["start"],
                    end=feature["end"],
                    score=feature["score"],
                    domain=feature["name"],
                    enzyme=enz_type,
                    ref=count,
                    locus=cds_id,
                    translation=feature["translation"]
                )
                count += 1

    def parse(self):
        count = 0
        for i, record in enumerate(gb_io.iter(self.args.i)):
            print(f"Reading Contig {i}", end='\r')
            self.process_record(record, count)
        return self.contigs

@jit(nopython=True)
def iupac_match_numba(char1, char2, iupac_lookup):
    return iupac_lookup[ord(char1), ord(char2)]

@jit(nopython=True)
def smith_waterman_jit(sequence1, sequence2, exact_match_score, match_score, mismatch_score, gap_penalty, iupac_lookup):
    # Initialize the score matrix and traceback matrix
    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    score_matrix = np.zeros((rows, cols))
    traceback_matrix = np.zeros((rows, cols), dtype=np.int32)

    # Fill the score and traceback matrices
    for i in range(1, rows):
        for j in range(1, cols):
            if sequence1[i - 1] == sequence2[j - 1]:  # If the characters are the same
                to_add = exact_match_score
            elif iupac_match_numba(sequence1[i - 1], sequence2[j - 1], iupac_lookup):  # If they are an IUPAC match
                to_add = match_score
            else:  # If they mismatch
                to_add = mismatch_score
            match = score_matrix[i - 1, j - 1] + to_add  # Match, so move up in both sequences
            delete = score_matrix[i - 1, j] + gap_penalty  # Move up one in sequence 2
            insert = score_matrix[i, j - 1] + gap_penalty  # Move up one in sequence 1
            score_matrix[i, j] = max(0, match, delete, insert)

            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = 1  # Diagonal
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = 2  # Up
            elif score_matrix[i, j] == insert:
                traceback_matrix[i, j] = 3  # Left

    # Find the cell with the maximum score in the score matrix
    max_score = np.max(score_matrix)
    max_i, max_j = 0, 0
    for i in range(rows):
        for j in range(cols):
            if score_matrix[i, j] == max_score:
                max_i, max_j = i, j
                break

    # Traceback to find the alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_i, max_j

    while i > 0 and j > 0 and score_matrix[i, j] > 0:
        if traceback_matrix[i, j] == 1:  # Diagonal
            aligned_seq1.append(sequence1[i - 1])
            aligned_seq2.append(sequence2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 2:  # Up, gap in sequence 2
            aligned_seq1.append(sequence1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i, j] == 3:  # Left, gap in sequence 1
            aligned_seq1.append('-')
            aligned_seq2.append(sequence2[j - 1])
            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return ''.join(aligned_seq1), ''.join(aligned_seq2), max_score

# Class definition remains mostly the same
class SmithWaterman:
    def __init__(self):
        # Define IUPAC codes
        self.iupac = {
            'A': 'A',
            'C': 'C',
            'G': 'G',
            'T': 'T',
            'R': 'AG',
            'Y': 'CT',
            'S': 'GC',
            'W': 'AT',
            'K': 'GT',
            'M': 'AC',
            'B': 'CGT',
            'D': 'AGT',
            'H': 'ACT',
            'V': 'ACG',
            'N': 'ACGT'
        }

        # Preprocess the IUPAC codes into a lookup table. 90 x 90 b/c ord(Y) is 89. Using ord makes indexing easier
        self.iupac_lookup = np.zeros((90, 90), dtype=np.bool_)

        for key, value in self.iupac.items():
            for char in value:
                self.iupac_lookup[ord(key), ord(char)] = True

    def smith_waterman(self, sequence1, sequence2, exact_match_score=3, match_score=1, mismatch_score=-1, gap_penalty=-2):
        aligned_seq1, aligned_seq2, max_score = smith_waterman_jit(sequence1, sequence2, exact_match_score, match_score, mismatch_score, gap_penalty, self.iupac_lookup)
        return aligned_seq1, aligned_seq2, max_score

def calculate_similarity(string1, string2):
    set1 = set(string1)
    set2 = set(string2)

    intersection = set1.intersection(set2)
    union = set1.union(set2)

    return len(intersection) / len(union)

@jit(nopython=True)
def distance(A_start, A_end, B_start, B_end, topology, contig_length):
    circular_distance = 0

    if topology == "circular":
        circular_distance = contig_length - (
                    max(A_end, B_end) - min(A_start, B_start))  # Total length - distance between start and end

    return min(max(0, B_start - A_end, A_start - B_end),
               circular_distance)  # Return max b/c one of the start/end differences will always be negative. Then if circ dist is less return that