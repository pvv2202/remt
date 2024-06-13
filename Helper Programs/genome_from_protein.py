from Bio import Entrez
from domainator.utils import parse_seqfiles
import subprocess
import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Annotated genbank file input')
parser.add_argument('-o', type=str, nargs='+', default=None, required=True, help='Output folder name')
parser.add_argument('--email', type=str, nargs='+', default=None, required=True, help='Email address for NCBI API access.')
args = parser.parse_args()

Entrez.email = args.email

undownloaded = []

def get_genome_accession_from_protein(protein_acc):
    try:
        # Link from protein to nucleotide
        handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_acc, linkname="protein_nuccore")
        records = Entrez.read(handle)
        handle.close()

        nuccore_ids = [link["Id"] for linksetdb in records for link in linksetdb["LinkSetDb"][0]["Link"]]
        if not nuccore_ids:
            return None

        # Link from nucleotide to assembly
        handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_ids, linkname="nuccore_assembly")
        records = Entrez.read(handle)
        handle.close()

        assembly_ids = [link["Id"] for linksetdb in records for link in linksetdb["LinkSetDb"][0]["Link"]]
        if not assembly_ids:
            return None

        # Fetch assembly summaries
        handle = Entrez.esummary(db="assembly", id=assembly_ids)
        records = Entrez.read(handle)
        handle.close()

        genome_accessions = []
        for record in records["DocumentSummarySet"]["DocumentSummary"]:
            if "LastMajorReleaseAccession" in record:
                genome_accessions.append(record["LastMajorReleaseAccession"])

        return genome_accessions if genome_accessions else None
    except Exception as e:
        print(f"An error occurred while processing {protein_acc}: {e}. Protein is likely TPA")
        undownloaded.append(protein_acc)
        return None

records = parse_seqfiles(args.i)

os.mkdir(args.o[0])
# for i, rec in enumerate(records):  # Each rec is a protein
#     protein_acc = rec.annotations["accessions"][0]
#     genome_acc = get_genome_accession_from_protein(protein_acc)
#     if genome_acc:
#         print(f"Protein Accession: {protein_acc}, Genome Accession: {', '.join(genome_acc)}")
#         os.mkdir(protein_acc)
#         for gen_ac in genome_acc:
#             command = f"datasets download genome accession {gen_ac}"
#             # Run the command
#             result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#
#             # Check if the command was successful
#             if result.returncode == 0:
#                 print("Download successful.")
#                 # Rename the downloaded file
#                 try:
#                     os.rename("ncbi_dataset.zip", f"{gen_ac}.zip")
#                     shutil.move(f"{gen_ac}.zip", protein_acc)
#                 except FileNotFoundError:
#                     print("Downloaded file not found. Make sure the download was successful.")
#                 except Exception as e:
#                     print(f"An error occurred while renaming the file: {e}")
#         try:
#             shutil.move(protein_acc, args.o[0])
#         except Exception as e:
#             print(f"An error occurred while moving the folder: {e}")
#     else:
#         print(f"Protein Accession: {protein_acc}, Genome Accession: Not found")

text = open(f"{args.o[0]}/undownloaded.txt", "w")
for item in undownloaded:
    text.write(item + "\n")
text.close()
