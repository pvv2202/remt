import argparse
import pickle
from domainator.seq_dist import seq_dist
from domainator.data_matrix import DataMatrix
from domainator.utils import parse_seqfiles, DomainatorCDS

parser = argparse.ArgumentParser(description='Generate HTML table from annotated GenBank file.')
parser.add_argument('-i', type=str, default=None, required=True, help='Annotated genbank or pickle file input')
parser.add_argument('-o', type=str, default="def", required=False, help='File output name (HTML)')
args = parser.parse_args()

# TODO: ADD ABILITY TO FEED IN GENERATED MATRICES
# TODO: MAKE EVERYTHING GO INTO A FOLDER
# TODO: PROBABLY MAKE A UTILS FILE TO START CLEANING SOME OF THIS UP
# TODO: FILTER SO THAT IT DOESNT TAKE AS LONG

class Contig:
    def __init__(self, ref, length, topology, strain):
        self.ref = ref
        self.length = length
        self.topology = topology
        self.strain = strain
        self.fusions = 0
        self.mts = {}
        self.res = {}
        self.hits = {}

class Enzyme:
    def __init__(self, seq, start, end, score, domain, enzyme, ref, locus, translation):
        # Reference number for faster lookup
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.domain = domain
        self.enzyme = enzyme
        self.ref = ref
        self.locus = locus
        self.translation = translation
        self.range = list()
        self.matches = list()

class Methyl(Enzyme):
    pass

class RE(Enzyme):
    pass

def distance(A_start, A_end, B_start, B_end, topology, contig_length):
    circular_distance = 0

    if topology == "circular":
        circular_distance = contig_length - (max(A_end, B_end) - min(A_start, B_start)) # Total length - distance between start and end

    return min(max(0, B_start - A_end, A_start - B_end), circular_distance) # Return max b/c one of the start/end differences will always be negative. Then if circ dist is less return that

prefixes = {"M.", "M1.", "M2.", "M3.", "M4."}

contigs = {}
re_sequences=""
mt_sequences=""

# Processing input file
if args.i.endswith('.gb'):
    # Populating contigs from input file
    print("Reading Input File")
    records = parse_seqfiles([args.i])

    print("Populating Contigs")
    for i, rec in enumerate(records):  # Each rec is a contig
        print(f"Reading Contig {i}", end='\r')
        count = 0
        cds_features = DomainatorCDS.list_from_contig(rec)

        for feature in cds_features:
            if len(feature.domain_features) == 0:
                continue

            domain_features = feature.domain_features[0]

            contig_id = rec.id
            if contig_id not in contigs:
                contigs[contig_id] = (
                    Contig(
                        ref=contig_id,
                        length=len(rec.seq),
                        topology=rec.annotations['topology'],
                        strain=rec.annotations['organism']
                    )
                )

            desc = domain_features.qualifiers['description'][0]
            score = domain_features.qualifiers['score'][0]

            # Filtering out to only consider applicable enzymes
            if "enzyme/methyltransferase" in desc:
                contigs[contig_id].fusions += 1
            elif all(keyword not in desc for keyword in ["control protein", "homing endonuclease", "subunit", "helicase", "nicking endonuclease", "orphan", "methyl-directed", "enzyme/methyltransferase"]) and "RecSeq:" in desc:
                # Extracting the RecSeq
                recseq_start = desc.find("RecSeq:") + len("RecSeq:")
                recseq_end = desc.find(";", recseq_start)
                seq = desc[recseq_start:recseq_end].split(", ")
                for s in seq:
                    if "-" in s:
                        seq.remove(s)
                        continue
                    if "NNN" in s:
                        seq.remove(s)

                if len(seq) == 0 or (len(seq) == 1 and "NNN" in seq[0]):
                    continue

                # Extracting the EnzType
                if "EnzType:" in desc and len(seq) > 0:
                    enz_start = desc.find("EnzType:") + len("EnzType:")
                    enz_end = desc.find(";", enz_start)
                    enz_type = desc[enz_start:enz_end]

                    # If MT/RE, respond accordingly - count is used as the reference number
                    name = domain_features.qualifiers['name'][0]
                    if any(name.startswith(prefix) for prefix in prefixes):
                        contigs[contig_id].mts[count] = (
                            Methyl(
                               seq=seq,
                               start=domain_features.location.start,
                               end=domain_features.location.end,
                               score=score,
                               domain=name,
                               enzyme=enz_type,
                               ref=count,
                               locus=domain_features.qualifiers['cds_id'][0],
                               translation=feature.feature.qualifiers['translation'][0]
                               )
                        )
                        mt_sequences += f">{domain_features.qualifiers['cds_id'][0]}\n{feature.feature.qualifiers['translation'][0]}\n"
                    else:
                        contigs[contig_id].res[count] = (
                            RE(
                               seq=seq,
                               start=domain_features.location.start,
                               end=domain_features.location.end,
                               score=score,
                               domain=name,
                               enzyme=enz_type,
                               ref=count,
                               locus=domain_features.qualifiers['cds_id'][0],
                               translation=feature.feature.qualifiers['translation'][0]
                               )
                        )
                        re_sequences += f">{domain_features.qualifiers['cds_id'][0]}\n{feature.feature.qualifiers['translation'][0]}\n"
                    count += 1

    # Save the contigs dictionary as a pickle object
    with open(f'{args.o}.pkl', 'wb') as f:
        pickle.dump(contigs, f)
elif args.i.endswith(".pkl"):
    with open(f'{args.i}', 'rb') as f:
        contigs = pickle.load(f)
        for contig in contigs.values():
            for r in contig.res.values():
                re_sequences += f">{r.locus}\n{r.translation}\n"
            for m in contig.mts.values():
                mt_sequences += f">{m.locus}\n{m.translation}\n"
else:
    print('Error: Unsupported file format. Please provide a .gb or .pkl file.')
    exit(1)

# Write fasta files, generate matrices, load matrices
with open(f"{args.o}_re_sequences.fasta", "w") as f:
    f.write(re_sequences)
with open(f"{args.o}_mt_sequences.fasta", "w") as f:
    f.write(mt_sequences)
seq_dist(f"{args.o}_re_sequences.fasta", "fasta", f"{args.o}_re_sequences.fasta", "fasta", None, "diamond_us", "score", 8, None, None, f"{args.o}_re_simlarity_matrix.hdf5", 0)
seq_dist(f"{args.o}_mt_sequences.fasta", "fasta", f"{args.o}_mt_sequences.fasta", "fasta", None, "diamond_us", "score", 8, None, None, f"{args.o}_mt_simlarity_matrix.hdf5", 0)
re_matrix = DataMatrix.from_file(f"{args.o}_re_simlarity_matrix.hdf5")
mt_matrix = DataMatrix.from_file(f"{args.o}_mt_simlarity_matrix.hdf5")

all_re = {r.locus: (c.ref, r) for c in contigs.values() for r in c.res.values()}
hits = {}

# Compare the systems
count = 0
for i, c in enumerate(contigs.values()):
    print(f'{int(i / len(contigs) * 100)}%', end='\r')
    for r in c.res.values():
        for m in c.mts.values():
            # If distance under threshold, check for matching re and then for matching mt within the right distance. If all applies, add to hits
            s1_dist = distance(min(r.start, r.end), max(r.start, r.end), min(m.start, m.end), max(m.start, m.end), c.topology, c.length)
            if s1_dist > 2500:
                continue
            for r2_loc in all_re:
                if re_matrix.data[re_matrix.row_to_idx[r.locus], re_matrix.column_to_idx[r2_loc]] < 50:
                    continue
                r2 = all_re[r2_loc][1]
                c2 = contigs[all_re[r2_loc][0]]
                for m2 in c2.mts.values(): # Check mts on the re's contig
                    if mt_matrix.data[mt_matrix.row_to_idx[m.locus], mt_matrix.column_to_idx[m2.locus]] < 50:
                        continue
                    s2_dist = distance(min(r2.start, r2.end), max(r2.start, r2.end), min(m2.start, m2.end), max(m2.start, m2.end), c2.topology, c2.length)
                    if s2_dist < 5000:
                        continue
                    hits[count] = {
                        "Contig": [c.ref, c2.ref],
                        "Strain": [c.strain, c2.strain],
                        "Types": [m.enzyme, r.enzyme, m2.enzyme, r2.enzyme],
                        "Domains": [m.domain, r.domain, m2.domain, r2.domain],
                        "Distance": [s1_dist, s2_dist],
                        "Sequences": [m.seq, r.seq, m2.seq, r2.seq],
                        "Loci": [m.locus, r.locus, m2.locus, r2.locus],
                        "Scores": [m.score, r.score, m2.score, r2.score],
                        "Translation": [m.translation, r.translation, m2.translation, r2.translation],
                    }
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
    html_table = "<table id='myTable' style='border-collapse: collapse; font-family: Courier New, monospace; font-size: 0.7em; width: 100%'>\n"
    # Create table header
    html_table += "<thead><tr>"
    for i, key in enumerate(header_keys):
        html_table += f"<th class='col-{i}' style='border: 1px solid black; padding: 8px; text-align: left;' onclick='toggleColumnWidth({i})'>" + str(key) + "</th>"
    html_table += "</tr></thead>\n"
    # Create table rows
    html_table += "<tbody>"
    for item in hits.values(): # 2 dictionaries
        if isinstance(item, dict):
            html_table += "<tr>"
            for i, (key, value) in enumerate(item.items()):
                cell_content = ""
                if key in ["Sequences", "Scores", "Types", "Domains", "Loci", "Translation"]:
                    cell_content = "R: " + str(value[1]) + "<br>" + "M: " + str(value[0]) + "<br>" + "R2: " + str(value[3]) + "<br>" + "M2: " + str(value[2])
                elif key in ["Contig", "Strain"]:
                    cell_content = "C1: " + str(value[0]) + "<br>" + "C2: " + str(value[1])
                else:
                    cell_content = "S1: " + str(value[0]) + "<br>" + "S2: " + str(value[1])
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