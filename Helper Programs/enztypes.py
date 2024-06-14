import argparse
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, default=None, required=True, help='File input (JSON or TSV)')
args = parser.parse_args()

#Handling input file
if args.i is not None:
    in_file = args.i
    if in_file.endswith('.json'):
        with open(in_file, 'r') as file:
            data = json.load(file)
    elif in_file.endswith('.tsv'):
        with open(in_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            data = [row for row in reader]
    else:
        print('Error: Unsupported file format. Please provide a JSON or TSV file.')
        exit(1)
else:
    print('Error: No input file specified.')
    exit(1)

enz_types = list()

for item in data:
    desc = item.get("domain_descriptions")
    if "EnzType:" in desc:
        i = desc.find("EnzType:")
        start = i + len("EnzType:")
        end = desc.find(";", start)
        type = desc[start:end]
    if type not in enz_types:
        enz_types.append(type)

print(enz_types)
