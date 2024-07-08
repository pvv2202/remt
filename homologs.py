'''
Sort through homologs of an MT and RE to output distances of systems and whether other enzymes are nearby
'''

import argparse
import os
from utils import Annotations, distance

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Restriction enzyme and methyltransferase directories')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
args = parser.parse_args()

class Homolog:
    '''
    Represents a homologous protein. Contains a dictionary of genomes in which it appears, each of which contains a dictionary of contig objects from utils.py.
    '''
    def __init__(self, id):
        self.genomes = {}
        self.id = id

class Genome:
    '''
    Represents a genome. Contains a dictionary of contig objects from utils.py.
    '''
    def __init__(self, id):
        self.contigs = {}
        self.contigs_all = {}
        self.systems = []
        self.id = id

re_homologs = {}
mt_homologs ={}

re_gs = []
mt_gs = []
re_con = []
mt_con = []

# RE directory, MT directory. Fill re_homologs and mt_homologs accordingly and find systems as you go
for i, dir in enumerate(args.i):
    for homolog_dir in os.listdir(dir):
        if not os.path.isdir(os.path.join(dir, homolog_dir)):
            continue
        dict = {}
        if i == 0 or i == 2:
            if homolog_dir not in re_homologs:
                re_homologs[homolog_dir] = Homolog(homolog_dir)
            dict = re_homologs
        elif i == 1 or i == 3:
            if homolog_dir not in mt_homologs:
                mt_homologs[homolog_dir] = Homolog(homolog_dir)
            dict = mt_homologs
        for genome in os.listdir(os.path.join(dir, homolog_dir)):
            if "annotated" not in genome:
                continue
            if i < 2:
                dict[homolog_dir].genomes[genome] = Genome(genome)
            # Reformat args to be passed to utils, then parse
            temp_args_dict = vars(args).copy()
            temp_args_dict['i'] = os.path.join(dir, homolog_dir, genome)
            temp_args = argparse.Namespace(**temp_args_dict)

            annotation = Annotations(temp_args)
            contigs = annotation.parse()

            for contig_id, contig in contigs.items():
                # Populate contigs_all with gs annotated contigs. Then find systems and reference back to see if anything is nearby
                if i < 2:
                    dict[homolog_dir].genomes[genome].contigs_all[contig_id] = contig
                    continue
                elif i >= 2:
                    dict[homolog_dir].genomes[genome].contigs[contig_id] = contig

                    # For each re, go through each mt and find the nearest one (for just the ones annotated with confirmed enzymes)
                    for re in contig.res.values():
                        nearest = [None, None, float('-inf'), None, None]
                        for mt in contig.mts.values():
                            dist = distance(re.start, re.end, mt.start, mt.end, contig.topology, contig.length)
                            if dist > nearest[2]:
                                nearest = [re, mt, dist, "None", "None"]
                        if nearest[0]:
                            if contig_id in dict[homolog_dir].genomes[genome].contigs_all:
                                mts_near = []
                                res_near = []
                                for mt2 in dict[homolog_dir].genomes[genome].contigs_all[contig_id].mts.values():
                                    if distance(nearest[0].start, nearest[0].end, mt2.start, mt2.end, contig.topology, contig.length) < 5000:
                                        mts_near.append(mt2.locus)
                                if len(mts_near) > 0:
                                    nearest[3] = mts_near
                                nearest.append(mts_near)
                                for re2 in dict[homolog_dir].genomes[genome].contigs_all[contig_id].res.values():
                                    if distance(nearest[1].start, nearest[1].end, re2.start, re2.end, contig.topology, contig.length) < 5000:
                                        res_near.append(re2.locus)
                                if len(res_near) > 0:
                                    nearest[4] = res_near
                            dict[homolog_dir].genomes[genome].systems.append(nearest)

dicts = [re_homologs, mt_homologs]
for dict in dicts:
    for homolog, homolog_obj in dict.items():
        for genome, genome_obj in homolog_obj.genomes.items():
            for system in genome_obj.systems:
                print(f"Homolog: {homolog}, Genome: {genome}, Distance: {system[2]}, Nearby MTs: {system[3]}, Nearby REs: {system[4]}")