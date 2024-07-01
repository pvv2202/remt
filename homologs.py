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
        self.systems = {}
        self.id = id

re_homologs = {}
mt_homologs ={}

# RE directory, MT directory. Fill re_homologs and mt_homologs accordingly and find systems as you go
for i, dir in enumerate(args.i):
    if not os.path.isdir(dir):
        continue
    for homolog_dir in os.listdir(dir):
        if not os.path.isdir(os.path.join(dir, homolog_dir)):
            continue
        dict = {}
        if i == 0:
            re_homologs[homolog_dir] = Homolog(homolog_dir)
            dict = re_homologs
        elif i == 1:
            mt_homologs[homolog_dir] = Homolog(homolog_dir)
            dict = mt_homologs
        for genome in os.listdir(os.path.join(dir, homolog_dir)):
            dict[homolog_dir].genomes[genome] = Genome(genome)
            # Reformat args to be passed to utils, then parse
            temp_args_dict = vars(args).copy()
            temp_args_dict['i'] = os.path.join(dir, homolog_dir, genome)
            temp_args = argparse.Namespace(**temp_args_dict)

            annotation = Annotations(temp_args)
            contigs = annotation.parse()

            for contig_id, contig in contigs.items():
                dict[homolog_dir].genomes[genome].contigs[contig_id] = contig
                count = 0
                # For each re, go through each mt and find the nearest one
                for re in contig.res.values():
                    nearest = (None, None, float('-inf'))
                    for mt in contig.mts.values():
                        dist = distance(re.start, re.end, mt.start, mt.end, contig.topology, contig.length)
                        if dist > nearest[2]:
                            nearest = (re, mt, dist)
                    if nearest[0]:
                        dict[homolog_dir].genomes[genome].systems[count] = nearest
                        count += 1

for homolog, homolog_obj in re_homologs.items():
    for genome, genome_obj in homolog_obj.genomes.items():
        for system, system_obj in genome_obj.systems.items():
            print(f"Homolog: {homolog}, Genome: {genome}, Distance: {system_obj[2]}")

for homolog, homolog_obj in mt_homologs.items():
    for genome, genome_obj in homolog_obj.genomes.items():
        for system, system_obj in genome_obj.systems.items():
            print(f"Homolog: {homolog}, Genome: {genome}, Distance: {system_obj[2]}")






