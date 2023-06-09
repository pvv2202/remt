import json
import argparse
import os
import csv

parser = argparse.ArgumentParser(description='Generate HTML table from JSON or TSV data.')
parser.add_argument('-i', type=str, default=None, required=True, help='File input (JSON or TSV)')
parser.add_argument('-o', type = str, default = "def", required = False, help='File output name (JSON and HTML)')
parser.add_argument('--coords', default = None, required = False, action='store_true', help='Include re/mt start and end points')
parser.add_argument('--shortest', default = None, required = False, action='store_true', help='Filter to only show the shortest re/mt pair for each exact match')
parser.add_argument('--shortest_all', default = None, required = False, action='store_true', help='Filter to show only the shortest re/mt pair for all matches on a contig for a specific re')
parser.add_argument('--ignore_short', default = None, required = False, action='store_true', help='Filter to ignore MTs with rec seq <3 bp')
parser.add_argument('--min_distance', default = None, required = False, type = int, help='Filter by distance')
parser.add_argument('--min_score', default = None, required = False, type = int, help='Filter by score')
parser.add_argument('--mt_in_range', default = None, required = False, type = int, help='Include a metric of how many mts are in a specified range of each re')
args = parser.parse_args()

#Handling input file
if args.i is not None:
    in_file = args.i
    if in_file.endswith('.json'):
        with open(in_file, 'r') as file:
            data = json.load(file)
    elif in_file.endswith('.tsv'):
        with open(in_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            data = [row for row in reader]
    else:
        print('Error: Unsupported file format. Please provide a JSON or TSV file.')
        exit(1)
else:
    print('Error: No input file specified.')
    exit(1)

class contig:
    def __init__(self, id):
        self.id = id
        self.mts = list()
        self.res = list()

class enzyme:
    def __init__(self, seq, start, end, score, domain, enzyme):
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.domain = domain
        self.enzyme = enzyme

class methyl(enzyme):
    def __init__(self, seq, start, end, score, domain, enzyme):
        super().__init__(seq, start, end, score, domain, enzyme)

class re(enzyme):
    def __init__(self, seq, start, end, score, domain, enzyme):
        super().__init__(seq, start, end, score, domain, enzyme)

def match(mseq, rseq):
    #Return false if the methyl sequence is greater than the re sequence
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

def filterAll(hits):
    updated_hits = list()
    to_remove = set()

    for i, check in enumerate(hits):
        for j, other in enumerate(hits):
            if i != j and check["Contig"] == other["Contig"] and check["RESeq"] == other["RESeq"]:
                #Adding the index of the hit to remove
                if check["Distance"] > other["Distance"]:
                    to_remove.add(i)
                else:
                    to_remove.add(j)
    
    updated_hits = [hit for i, hit in enumerate(hits) if i not in to_remove]

    return updated_hits

def distance(A_start, A_end, B_start, B_end):
    return max(0, B_start - A_end, A_start - B_end)

contigs = list()

prefixes = {"M.", "M1.", "M2.", "M3.", "M4."}

for item in data:
    isRE = True
    not_in_contigs = True
    desc = item.get("domain_descriptions")
    if all(keyword not in desc for keyword in ["control", "homing", "enzyme/methyltransferase", "subunit", "helicase", "nicking"]) and "RecSeq:" in desc:        #Add contigs
        #Determining which contig to add to. If not created, create it and append to list
        for c in contigs:
            if item.get("contig") == c.id:
                ctg = c
                not_in_contigs = False
        if not_in_contigs:
            ctg = contig(item.get("contig"))
            contigs.append(ctg)
        i = desc.find("RecSeq:")
        recseq_start = i + len("RecSeq:")
        recseq_end = desc.find(";", recseq_start)
        seq = desc[recseq_start:recseq_end]
        if not args.ignore_short or len(seq) > 2:
            if "EnzType:" in desc:
                j = desc.find("EnzType:")
                et_start = j + len("EnzType:")
                et_end = desc.find(";", et_start)
                et = desc[et_start:et_end]
                for m in prefixes:
                    if item.get("domain").startswith(m):
                        ctg.mts.append(methyl(seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain"), et))
                        isRE = False
            if isRE:
                ctg.res.append(re(seq, item.get("start"), item.get("end"), item.get("score"), item.get("domain"), et))

file.close()

hits = list()

for c in contigs:
    for m in c.mts:
        for r in c.res:
            if match(m.seq, r.seq):
                mt_count = 0
                if args.mt_in_range is not None:
                    mt_count = sum(1 for mt_in_range in c.mts if distance(min(r.start, r.end), max(r.start, r.end), min(mt_in_range.start, mt_in_range.end), max(mt_in_range.start, mt_in_range.end)) < args.mt_in_range)
                temp = {
                    "Contig": c.id,
                    "REType": r.enzyme,
                    "MT": m.domain,
                    "RE": r.domain,
                    "MTSeq": m.seq,
                    "RESeq": r.seq,
                    "AvgScore": round((r.score + m.score) / 2),
                    "MTScore": m.score,
                    "REScore": r.score,
                    "Distance": distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end))
                }
                if args.mt_in_range is not None:
                    temp["MT within " + str(args.mt_in_range) + " bp"] = mt_count

                if args.coords:
                    temp["M S/E"] = str(m.start) + ", " + str(m.end)
                    temp["R S/E"] = str(r.start) + ", " + str(r.end)

                hits.append(temp)

#Argument Handling
if args.shortest:
    hits = filterShortest(hits)

if args.shortest_all:
    hits = filterAll(hits)

if args.min_distance is not None:
   hits = [h for h in hits if h["Distance"] >= args.min_distance]

if args.min_score is not None:
   hits = [h for h in hits if h["AvgScore"] >= args.min_score]

hit_data = json.dumps(hits, indent=4)

# Write the JSON hit data to a file
json_filename = args.o + '.json'
with open(json_filename, 'w') as file:
    file.write(hit_data)

# Check if the HTML file already exists
html_filename = args.o + ".html"
if os.path.exists(html_filename):
    os.remove(html_filename)

# Read the JSON data from the file
with open(json_filename, 'r') as file:
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

    # Write the HTML table to the file
    with open(html_filename, "a") as file:
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
    with open(html_filename, "a") as file:
        file.write(html_script)
else:
    print("Error: JSON File Empty")
