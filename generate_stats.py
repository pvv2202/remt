'''
Generate stats from remt2 HTML output. In newer versions or remt, this doesn't tell much because it only outputs split systems.
REMT can be easily modified to output full systems, however.
'''

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser(description='Generate statistics from remt2 HTML output.')
parser.add_argument('-i', type=str, nargs='+', default=None, required=True, help='Annotated genbank file input')
args = parser.parse_args()

if args.i is not None:
    with open(args.i[0], 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')

    table = soup.find('table', {'id': 'myTable'})

    headers = [header.text for header in table.find('thead').find_all('th')]
    distance_index = headers.index('Distance')

    distances = []

    for row in table.find('tbody').find_all('tr'):
        cells = row.find_all('td')
        distance_value = cells[distance_index].get_text(strip=True)
        distances.append(int(distance_value))
else:
    print('Error: No input file specified.')
    exit(1)

# Plot histogram
plt.hist(distances, bins=10, alpha=0.75, edgecolor='black', density=False)
plt.title('System Distances')
plt.xlabel('Distance')
plt.ylabel('Count')
plt.show()

# Plot KDE
sns.kdeplot(distances, fill=True, clip=(0, None))
plt.title('KDE of Distances')
plt.xlabel('Distance')
plt.ylabel('Density')
plt.show()