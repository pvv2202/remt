'''
Program to clean the output of compare systems to be used as an input for REMT
'''

import csv
import math

input_file = "compgen_comp.csv"
output_file = "cg_cdd.csv"

# Function to check if a value is NaN
def is_nan(value):
    try:
        return math.isnan(float(value))
    except:
        return False

# Read the original file and write only the selected columns to a new file
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile)
    fieldnames = ['Contig', 'Domains', 'Distance']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in reader:
        # Handle multi-line entries by cleaning up the fields
        contig = row['Contig']
        domains = row['Domains']
        distance = row['Distance']
        
        # Remove rows where "Domains", "Distance", or "Contig" is None or NaN
        if domains and not is_nan(domains):
            domains = domains.replace('\n', ' ').strip()
        else:
            continue
        
        if distance and not is_nan(distance):
            distance = distance.replace('\n', ' ').strip()
        else:
            continue

        if contig and not is_nan(contig):
            contig = contig.replace('\n', ' ').strip()
        else:
            continue

        writer.writerow({'Contig': contig, 'Domains': domains, 'Distance': distance})

print(f"Selected columns written to {output_file}")
