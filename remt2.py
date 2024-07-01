import argparse
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from utils import Annotations, SmithWaterman, distance, calculate_similarity

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Annotated genbank or pickle file input')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
parser.add_argument('--stats', default=None, required=False, action='store_true', help='Print stats to console')
parser.add_argument('--histogram', nargs='?', default=None, const=10, type=int, help='Output histogram of system distances. Specify bin size after. Default is 10')
parser.add_argument('--kde', default=None, required=False, action='store_true', help='Output kde of system distances (clipped to display positive values only)')
parser.add_argument('--excel', default=None, required=False, action='store_true', help='Make results excel friendly (default is HTML friendly)')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp or many Ns in rec seq (nonspecific)')
parser.add_argument('--min_distance', default=None, required=False, type=int, help='Filter by distance')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--in_range', default=5000, required=False, type=int, help='Include a metric of how many mts/res are in range of each hit')
args = parser.parse_args()

contigs = {}
sw = SmithWaterman()

if args.i.endswith(".gb"):
    annotations = Annotations(args)
    contigs = annotations.parse()

    # Save the contigs dictionary as a pickle object
    with open(f'{args.o}.pkl', 'wb') as f:
        pickle.dump(contigs, f)
elif args.i.endswith(".pkl"):
    with open(f'{args.i}', 'rb') as f:
        contigs = pickle.load(f)
else:
    print('Error: Unsupported file format. Please provide a .gb or .pkl file.')
    exit(1)

print("\r")
print("Finding Systems")
for i, c in enumerate(contigs.values()):
    print(f'{int(i / len(contigs) * 100)}%', end='\r')
    # For each restriction enzyme and methyltransferase
    for r in c.res.values():
        for rseq in r.seq:
            for m in c.mts.values():
                for mseq in m.seq:
                    dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length)
                    # Add any re/mts within specified range
                    if args.in_range and dist < args.in_range:
                        m.range.append(f'{r.domain}: {r.seq}')
                        r.range.append(f'{m.domain}: {m.seq}')

                    # If mseq <= rseq (if greater cannot block every site), run smith-waterman
                    if len(mseq) <= len(rseq):
                        # Filter any weird matches to nonsepc regions
                        if ("NNN" in rseq and "NNN" not in mseq) or ("NNN" in mseq and "NNN" not in rseq):
                            continue
                        aligned_mseq, aligned_rseq, score = sw.smith_waterman(mseq, rseq)
                        #Filtering by score, length of rseq, absolute value of mseq
                        if score > 2 and len(aligned_rseq) > 3 and (len(aligned_mseq) == len(mseq)):
                            sim = calculate_similarity(mseq, rseq)
                            asim = calculate_similarity(aligned_mseq, aligned_rseq)

                            #Doing a better job getting the score
                            weighted_score = round((asim + sim) * score) + round(1000/dist) if dist > 0 else 1000

                            if r.ref not in c.hits or weighted_score > c.hits[r.ref]["Score"]:
                                temp = {
                                    "Contig": c.ref,
                                    "Strain": c.strain,
                                    "Types": [m.enzyme, r.enzyme],
                                    "Domains": [m.domain, r.domain],
                                    "Sequences": [mseq, rseq],
                                    "Alignment": [aligned_mseq, aligned_rseq],
                                    "Loci": [m.locus, r.locus],
                                    "Score": weighted_score,
                                    "Scores": [m.score, r.score],
                                    "Distance": dist,
                                    "Coords": [str(m.start) + ", " + str(m.end), str(r.start) + ", " + str(r.end)],
                                    "Translation": [m.translation, r.translation]
                                }

                                if args.in_range:
                                    if len(r.range) == 0:
                                        r.range.append("None")
                                    if len(m.range) == 0:
                                        m.range.append("None")
                                    for mt_in_range in r.range:
                                        temp["MT within " + str(args.in_range) + " bp"] = mt_in_range
                                    for rt_in_range in m.range:
                                        temp["RE within " + str(args.in_range) + " bp"] = rt_in_range

                                c.hits[r.ref] = temp

print("\r    ")

hits = {}
count = 0
for contig in contigs.values():
    for hit in contig.hits.values():
        hits[count] = hit
        count += 1

#Filtering hits for distance
if args.min_distance is not None:
    hits = {h: v for h, v in hits.items() if v["Distance"] >= args.min_distance}

print(len(hits))

#Printing necessary data for stats
if args.stats:
    num_none = 0
    num_both = 0
    num_only_re = 0
    num_only_mt = 0
    hit_stats = []
    re_only = []

    for c in contigs.values():
        if (len(c.res) > 0 and len(c.mts) > 0) or c.fusions > 0:
            num_both += 1
            curr = c.fusions + len(c.hits)
            if len(hit_stats) < curr + 1:
                hit_stats.extend([0] * (curr - len(hit_stats) + 1))
            hit_stats[curr] += 1
        elif len(c.res) == 0 and len(c.mts) == 0:
            num_none += 1
        elif len(c.res) > 0 and len(c.mts) == 0:
            num_only_re += 1
            if c not in re_only:
                re_only.append(c)
        elif len(c.res) == 0 and len(c.mts) > 0:
            num_only_mt += 1

    print("Data for Contigs with both MT and RE:")
    for i, h in enumerate(hit_stats):
        print(f"Contigs with {i} Systems: {h}")

    print("\nData for Contigs in general:")
    print("Num Contigs: " + str(len(contigs)))
    print("Num Both/Fusion: " + str(num_both))
    print("Num None: " + str(num_none))
    print("Num Only RE: " + str(num_only_re))
    print("Num Only MT: " + str(num_only_mt))

if args.histogram is not None or args.kde:
    distances = [hit["Distance"] for hit in hits.values()]
    if args.histogram:
        # Plot histogram
        plt.hist(distances, bins=args.histogram, alpha=0.75, edgecolor='black', density=False)
        plt.title('System Distances')
        plt.xlabel('Distance')
        plt.ylabel('Count')
        plt.show()
    if args.kde:
        # Plot KDE
        sns.kdeplot(distances, fill=True, clip=(0, None))
        plt.title('KDE of Distances')
        plt.xlabel('Distance')
        plt.ylabel('Density')
        plt.show()

if hits:
    # Get the keys from the first dictionary in hits
    for key, value in hits.items():
        if isinstance(value, dict):
            header_keys = value.keys()
            break
    else:
        # Handle the case when no dictionaries are found in hits
        header_keys = []

    # Generate the HTML table with centered content
    html_table = "<table id='myTable' style='border-collapse: collapse; font-family: Courier New, monospace; font-size: 0.7em; width: 100%'>\n"
    # Create table header
    html_table += "<thead><tr>"
    for i, key in enumerate(header_keys):
        html_table += f"<th class='col-{i}' style='border: 1px solid black; padding: 8px; text-align: left;' onclick='toggleColumnWidth({i})'>" + str(key) + "</th>"
    html_table += "</tr></thead>\n"
    # Create table rows
    html_table += "<tbody>"
    for item in hits.values():
        if isinstance(item, dict):
            html_table += "<tr>"
            for i, (key, value) in enumerate(item.items()):
                cell_content = ""
                if key in ["Alignment", "Sequences", "Coords", "Scores", "Types", "Domains", "Loci", "Translation"]:
                    if args.excel:
                        cell_content = "R: " + str(value[1]) + "&#10;" + "M: " + str(value[0])
                    else:
                        cell_content = "R: " + str(value[1]) + "<br>" + "M: " + str(value[0])
                else:
                    cell_content = str(value)
                html_table += f"<td class='col-{i}' style='border: 1px solid black; padding: 8px; text-align: left;'>{cell_content}</td>"
            html_table += "</tr>\n"
    html_table += "</tbody></table>"

    output_filename = args.o + ".html"

    # Write the HTML table to the file
    with open(output_filename, 'w') as file:
        file.write(html_table)

    html_script = """
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.3/css/jquery.dataTables.min.css">
    <script src="https://cdn.datatables.net/1.11.3/js/jquery.dataTables.min.js"></script>
    <style>
        .collapsed-column {
            width: 10px !important;
            max-width: 10px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }
        table.dataTable td.collapsed-column, table.dataTable th.collapsed-column {
            width: 10px;
            max-width: 10px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }
    </style>
    <script>
        $(document).ready(function() {
            $('#myTable').DataTable();

            // Add click event to the whole table to toggle column width
            $('#myTable').on('click', 'td, th', function() {
                var index = $(this).index();
                toggleColumnWidth(index);
            });
        });

        function toggleColumnWidth(index) {
            var columnClass = '.col-' + index;
            $(columnClass).toggleClass('collapsed-column');
        }
    </script>
    """
    # Append the JavaScript code to the HTML file
    with open(output_filename, "a") as file:
        file.write(html_script)
