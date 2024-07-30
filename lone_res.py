import argparse
import pickle
from utils import Annotations, distance

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Annotated genbank or pickle file input')
parser.add_argument('--min_score', default=None, required=False, type=int, help='Filter by score (for enzymes, not for final weighted score)')
parser.add_argument('--ignore_nonspec', default=None, required=False, action='store_true', help='Filter to ignore MTs with rec seq < 3 bp (nonspecific)')
parser.add_argument('-o', type=str, default=None, required=True, help='File output name (HTML)')
args = parser.parse_args()

contigs = {}

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

    # for mt in c.mts.values():
    #     print(f"MT: {mt.enzyme}, Start: {mt.start}, End: {mt.end}")
    #
    # for re in c.res.values():
    #     print(f"RE: {re.enzyme}, Start: {re.start}, End: {re.end}")

    # For each restriction enzyme and methyltransferase
    if len(c.mts) > 50:
        continue

    for r in c.res.values():
        if "Type III" in r.enzyme:
            continue

        found_close = False
        for m in c.mts.values():
            dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length)

            if dist <= 10000:
                found_close = True
                break

        if found_close:
            continue

        c.hits[r.ref] = {
            "Contig": c.ref,
            "Strain": c.strain,
            "Type": r.enzyme,
            "Domain": r.domain,
            "Sequence": r.seq,
            "Loci": r.locus,
            "Scores": r.score,
            "Coords": str(r.start) + ", " + str(r.end),
            "Translation": r.translation
        }

print("\r    ")

hits = {}
count = 0
for contig in contigs.values():
    for hit in contig.hits.values():
        hits[count] = hit
        count += 1

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
    html_table = """
    <table id='myTable' style='border-collapse: collapse; font-family: Courier New, monospace; font-size: 0.7em; width: 100%'>
    <thead><tr>
    """
    # Create table header
    for i, key in enumerate(header_keys):
        html_table += f"<th class='col-{i}' style='border: 1px solid black; padding: 8px; text-align: left;' onclick='toggleColumnWidth({i})'>" + str(key) + "</th>"
    html_table += "</tr></thead>\n"
    # Create table rows
    html_table += "<tbody>"
    for item in hits.values():
        if isinstance(item, dict):
            html_table += "<tr>"
            for i, (key, value) in enumerate(item.items()):
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
            $('#myTable').DataTable({
                "paging": true,
                "searching": true,
                "ordering": true,
                "columnDefs": [{
                    "targets": "_all",
                    "orderable": true
                }]
            });
        });

        function toggleColumnWidth(columnIndex) {
            var column = $('#myTable').DataTable().column(columnIndex);
            $(column.header()).toggleClass('collapsed-column');
            column.nodes().to$().toggleClass('collapsed-column');
        }
    </script>
    """
    # Append the JavaScript code to the HTML file
    with open(output_filename, "a") as file:
        file.write(html_script)