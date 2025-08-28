from Bio import SeqIO

# Input files
file1 = "Csa_annotated_results_dedup.gb"
file2 = "Csa_hit_annotated.gb"
output_file = "Csa_annotated_results_dedup_merged.gb"

# Read records from both files
records1 = list(SeqIO.parse(file1, "genbank"))
records2 = list(SeqIO.parse(file2, "genbank"))

# Merge them into one list
merged_records = records1 + records2

# Write to a new GenBank file
SeqIO.write(merged_records, output_file, "genbank")

print(f"Merged {len(records1)} + {len(records2)} records into {output_file}")
