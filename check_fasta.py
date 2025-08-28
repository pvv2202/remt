import argparse

parser = argparse.ArgumentParser(description='Check 2 fasta files to see if there are any errors')
parser.add_argument('-check', type=str, default=None, required=True, help='Fasta file to check')
parser.add_argument('-ref', type=str, default=None, required=True, help='Reference fasta file to check against')
args = parser.parse_args()

def read_fasta_names(filepath):
    names = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Remove the ">" and strip whitespace
                name = line[1:].strip()
                names[name[3:]] = name[:3] # Name to characterization

    return names

# Example usage
check = read_fasta_names(args.check)
ref = read_fasta_names(args.ref)

for name in check.keys():
    if name not in ref:
        print(f"{name} not included!")
        continue

    if check[name] != ref[name]:
        print(f"{name} is currently {check[name]} but should be {ref[name]}")


