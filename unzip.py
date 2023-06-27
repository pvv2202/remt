import os
import gzip
import shutil

# specify directory where the files are
directory = 'to_zip'

# the name of the output file
output_file_name = 'concatenated.gbff'

with open(output_file_name, 'wb') as wfd:
    for filename in os.listdir(directory):
        if filename.endswith('.gbff.gz'):
            file_path = os.path.join(directory, filename)
            with gzip.open(file_path, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)