def convert_to_fasta(input_file, output_file):
    with open(input_file, 'r') as f_input, open(output_file, 'w') as f_output:
        line_count = 0
        for line in f_input:
            line = line.strip()
            if line:
                if line_count == 0:
                    f_output.write(f">{line}\n")
                else:
                    f_output.write(f"{line}\n")
                line_count += 1

# Specify the paths to the input text file and output FASTA file
input_file_path = 'convert.txt'
output_file_path = 'BsuRI.fasta'

# Call the conversion function
convert_to_fasta(input_file_path, output_file_path)