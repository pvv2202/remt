import pandas as pd

data = pd.read_csv("bacillus_contig.tsv", sep='\t')

count = 0
for item in data.iterrows():
    count += 1

print(count)

