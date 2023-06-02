import json
import argparse

parser = argparse.ArgumentParser(description='Generate HTML table from JSON data.')
parser.add_argument('-i', type = str, default = None, required = True, help='File input (JSON)')
parser.add_argument('-o', type = str, default = "def", required = False, help='File output name (JSON and HTML)')
parser.add_argument('--coords', default = None, required = False, action='store_true', help='Include re/mt start and end points')
parser.add_argument('--shortest', default = None, required = False, action='store_true', help='Filter to only show the shortest re/mt pair for each exact match')
parser.add_argument('--shortest_all', default = None, required = False, action='store_true', help='Filter to show only the shortest re/mt pair for all matches on a contig for a specific re')
parser.add_argument('--min_distance', default = None, required = False, type = int, help='Filter by distance')
parser.add_argument('--min_score', default = None, required = False, type = int, help='Filter by score')
args = parser.parse_args()

#Handling input file
if args.i is not None:
    in_file = args.i
    with open(in_file, 'r') as file:
        data = json.load(file)
else:
    print(FileExistsError)

class enzyme:
    def __init__(self, contig, seq, start, end, score, domain):
        self.contig = contig
        self.seq = seq
        self.start = start
        self.end = end
        self.score = score
        self.domain = domain

    def __str__(self):
        return f"Contig: {self.contig}, Seq: {self.seq}, Start: {self.start}, End: {self.end}, Score: {self.score}, Domain: {self.domain}"

class methyl(enzyme):
    def __init__(self, contig, seq, start, end, score, domain):
        super().__init__(contig, seq, start, end, score, domain)

class re(enzyme):
    def __init__(self, contig, seq, start, end, score, domain):
        super().__init__(contig, seq, start, end, score, domain)

def match(mseq, rseq):
    # Return false if the methyl sequence is greater than the re sequence
    if len(mseq) > len(rseq):
        return False
    else: 
        code = {
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['T'],
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T']
        }
        for i in range(len(rseq) - len(mseq) + 1):
            count = 0
            for j, char in enumerate(mseq.upper()):
                if char in code:
                    if rseq.upper()[i+j] in code[char]:
                        count += 1
            if count == len(mseq):
                return True

        return False
    
#Filter so that only the shortest distance re/methyl match is provided
def filterShortest(hits):
    updated_hits = list()
    to_remove = set()

    for i, check in enumerate(hits):
        for j, other in enumerate(hits):
            if i != j and check["Contig"] == other["Contig"] and check["MethylSeq"] == other["MethylSeq"] and check["RESeq"] == other["RESeq"]:
                #Adding the index of the hit to remove
                if check["Distance"] > other["Distance"]:
                    to_remove.add(i)
                else:
                    to_remove.add(j)
    
    updated_hits = [hit for i, hit in enumerate(hits) if i not in to_remove]

    return updated_hits

def expand_iupac(seq):
    code = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

    sequences = ['']
    for c in seq.upper():
        if c in code:
            sequences = [s + n for s in sequences for n in code[c]]
        else:
            sequences = [s + c for s in sequences]

    return sequences

methyls = list()
res = list()

for item in data:
    desc = item.get("domain_descriptions")
    if "RecSeq:" in desc:
        seq_list = list()
        i = desc.find("RecSeq:")
        recseq_start = i + len("RecSeq:")
        recseq_end = desc.find(";", recseq_start)
        seq = desc[recseq_start:recseq_end]
        if item.get("domain")[0] == 'M':
            methyls.append(methyl(item.get("contig"), seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain")))
        else:
            res.append(re(item.get("contig"), seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain")))

file.close()

hits = list()

#TODO: ADD ABILITY TO FILTER FOR MT/RE GROUPS WITHIN CONTIGS IN GENERAL SO THAT METHYLS WITH DIFFERENT SEQUENCES THAT BOTH BLOCK AN RE ARE REMOVED

#Loop through each methyl, if they have same rec seq and on same contig, add to hits. 
for m in methyls:
    for r in res:
        if match(m.seq, r.seq) and r.contig == m.contig:
            if args.coords:
                temp = {
                    "Contig": m.contig, "Methyl": m.domain, "RE": r.domain, "M S/E": str(m.start) + ", " + str(m.end), "R S/E": str(r.start) + ", " + str(r.end), "Distance": min(abs(r.end - m.start), abs(m.end - r.start)), "MethylSeq":m.seq, "RESeq":r.seq, "AvgScore":round((r.score + m.score)/2), 
                }
            else:
                temp = {
                    "Contig": m.contig, "Methyl": m.domain, "RE": r.domain, "Distance": min(abs(r.end - m.start), abs(m.end - r.start)), "MethylSeq":m.seq, "RESeq":r.seq, "AvgScore":round((r.score + m.score)/2), 
                }
            hits.append(temp)

#Argument handling
if args.shortest:
    hits = filterShortest(hits)

if args.min_distance is not None:
   hits = [h for h in hits if h["Distance"] >= args.min_distance]

if args.min_score is not None:
   hits = [h for h in hits if h["AvgScore"] >= args.min_score]

hit_data = json.dumps(hits, indent=4)

with open(args.o + '.json', 'w') as file:
    file.write(hit_data)

# Read the JSON data from the file
with open(args.o + '.json', 'r') as file:
    json_data = json.load(file)

if hit_data:
    # Generate the HTML table with centered content
    html_table = "<table id='myTable' style='border-collapse: collapse;'>\n"
    # Create table header
    html_table += "<thead><tr>"
    for key in json_data[0].keys():
        html_table += "<th style='border: 1px solid black; padding: 8px; text-align: left;'>" + key + "</th>"
    html_table += "</tr></thead>\n"
    # Create table rows
    html_table += "<tbody>"
    for item in json_data:
        html_table += "<tr>"
        for value in item.values():
            html_table += "<td style='border: 1px solid black; padding: 8px; text-align: left;'>" + str(value) + "</td>"
        html_table += "</tr>\n"
    html_table += "</tbody></table>"

    # Write the HTML table to a file
    with open(args.o + ".html", "a") as file:
        file.write(html_table)

    # Append the script tags to include the jQuery and DataTables libraries
    html_script = """
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.3/css/jquery.dataTables.min.css">
    <script src="https://cdn.datatables.net/1.11.3/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready(function() {
            $('#myTable').DataTable();
        });
    </script>
    """
    # Append the JavaScript code to the HTML file
    with open(args.o + ".html", "a") as file:
        file.write(html_script)
else:
    print("Error: JSON File Empty")
