import os
import subprocess
import glob

# Define the root directory where the protein directories are located
root_dir = "/Users/pvandervort/Downloads/BpeH640ORF18650P"
reference_fasta = "/Users/pvandervort/Downloads/Gold_Standard_np_prot.fasta"

# Function to annotate a GenBank file
def annotate_genbank(gb_file, reference_fasta, max_overlap=0.6):
    output_file = gb_file.replace(".gb", "_annotated.gb")
    command = [
        "conda", "run", "-n", "domainator", "domainate.py",
        "-i", gb_file,
        "-r", reference_fasta,
        "--max_overlap", str(max_overlap),
        "-o", output_file
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Annotated file saved as: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while annotating {gb_file}: {e}")

# Iterate over each protein directory and process the GenBank files
for i, protein_dir in enumerate(os.listdir(root_dir)):
    full_path = os.path.join(root_dir, protein_dir)
    if os.path.isdir(full_path):
        # Find all .gb files in the current protein directory
        gb_files = glob.glob(os.path.join(full_path, "*.gb"))
        for gb_file in gb_files:
            if "_annotated" not in gb_file:
                if gb_file.replace(".gb", "_annotated.gb") not in gb_files:
                    annotate_genbank(gb_file, reference_fasta)
                else:
                    print(f"File {gb_file} has already been annotated.")


