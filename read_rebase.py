'''
Program to read REBASE data
'''
import pandas as pd

class Contig:
    def __init__(self, ref):
        self.ref = ref
        self.mts = {}
        self.res = {}

class Enzyme:
    def __init__(self, seq, start, end):
        self.seq = set([seq])
        self.start = int(start)
        self.end = int(end)

class Methyl(Enzyme):
    pass

class RE(Enzyme):
    pass

def distance(start1, end1, start2, end2):
    return min(abs(start1 - end2), abs(start2 - end1))

def main():
    file_path = '~/Downloads/enzyme_info_for_genbank_accessions.2024-08-01.txt'
    prefixes = ('M.', 'M1.', 'M2.', 'M3.', 'M4.')

    # Read the data into a DataFrame
    data = pd.read_csv(file_path, sep='\t', header=None)
    data.columns = ['Num1', 'Name', 'Type', 'Seq', 'Num2', 'Name2', 'Contig', 'Start', 'End']
    data = data[['Name', 'Seq', 'Contig', 'Start', 'End']]

    contigs = {}

    # Populate contigs
    for index, enzyme in data.iterrows():
        name, seq, contig_id, start, end = enzyme['Name'], enzyme['Seq'], enzyme['Contig'], enzyme['Start'], enzyme['End']
        if contig_id in contigs:
            contig = contigs[contig_id]
        else:
            contig = Contig(contig_id)
            contigs[contig_id] = contig

        if '.' in name:
            if name.startswith(prefixes):
                if name not in contig.mts:
                    contig.mts[name] = Methyl(seq, start, end)
                else:
                    contig.mts[name].seq.add(seq)
            else:
                continue
        else:
            if name not in contig.res:
                contig.res[name] = RE(seq, start, end)
            else:
                contig.res[name].seq.add(seq)

    hits = {}
    count = 0
    for contig in contigs.values():
        for rid, r in contig.res.items():
            best = [r, None, float('-inf')]
            for mtid, m in contig.mts.items():
                dist = distance(r.start, r.end, m.start, m.end)
                if dist < 10000:
                    best[1] = None # Found a nearby MT, so ignore
                    break
                elif set(r.seq) & set(m.seq):
                    best = [r, m, dist]
            if best[1] != None:
                hits[count] = {
                    "Contig": contig.ref,
                    "RE": rid,
                    "MT": mtid,
                    "RSeq": best[0].seq,
                    "MSeq": best[1].seq,
                    "Distance": best[2],
                }
                count += 1

    for hit in hits.values():
        print(hit)
        print('\n')

if __name__ == "__main__":
    main()
