'''
Program to standardize the REBASE protein sequences to the Gold Standard protein sequences. Highly specialized for np_prot but could easily be modified
for a different FASTA file format.
'''

from Bio import SeqIO

to_change = "reb_all_np_prot.faa"
to_check = "All_REBASE_Gold_Standards_Protein_20221213.fasta"

mp = {"M.", "M1.", "M2.", "M3.", "M4.", "M5.", "M6."}
sp = {"S.", "S1.", "S2.", "S3.", "S4.", "S5."}

reb1 = SeqIO.parse(to_change, "fasta")  #rebase
reb2 = SeqIO.parse(to_change, "fasta")  #rebase
gs = SeqIO.parse(to_check, "fasta")  #gold standard

ds1 = {}
ds2 = {}

count = 0
repeats = 0

class Enz:
    def __init__(self, id, rseq):
        self.id = id
        self.rseq = []
        self.rseq.append(rseq)

for sequence in reb1:
    description_parts = sequence.description.split(' ', 1)
    description = description_parts[1]
    if sequence.seq in ds1.keys():
        repeats += 1
        if len(ds1[sequence.seq].rseq) == 1 and ds1[sequence.seq].rseq[0] == "-":
            ds1[sequence.seq].rseq[0] = description
        elif all(description != d for d in ds1[sequence.seq].rseq):
            ds1[sequence.seq].rseq.append(description)
    else:
        ds1[sequence.seq] = Enz(sequence.id, description)

print("Repeats: " + str(repeats))

for sequence in gs:
    ds2[sequence.seq] = sequence.id

seq2_records = []
iterated = {}

for sequence in reb2:
    description_parts = sequence.description.split(' ', 1)
    desc = description_parts[1]
    seq = sequence.seq
    id = sequence.id
    if seq in ds1.keys():
        if len(ds1[seq].rseq) == 1:
            desc = ds1[seq].rseq[0]
        else:
            for i, d in enumerate(ds1[seq].rseq):
                desc += d
                if i < len(ds1[seq].rseq) - 1:
                    desc += ", "
        id = ds1[seq].id

    if seq not in iterated.keys():
        if seq not in ds2.keys():
            iterated[seq] = True
            sequence.seq = seq.replace(" ", "")

            if any(id.startswith(prefix) for prefix in mp):
                toAdd = "methyltransferase"
            elif any(id.startswith(prefix) for prefix in sp):
                toAdd = "specificity subunit"
            elif id.startswith("C."):
                toAdd = "control protein"
            else:
                toAdd = "-"

            sequence.description = "EnzType:" + toAdd + "; " + "RecSeq:" + desc + "; " + "SeqLength:" + str(len(sequence.seq))

            seq2_records.append(sequence)
        else:
            count += 1

print("Overlap: " + str(count))

output_file = "reb_gs.fasta"
SeqIO.write(seq2_records, output_file, "fasta")
