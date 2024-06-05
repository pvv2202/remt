# remt

Program to find restriction modification systems. Optimized for finding high-confidence split systems.
This repository contains both remt.py and remt2.py. remt.py was a first pass at making this program and is potentially useful as a reference or for 
very specific needs, but remt2 is faster and more effective.

remt2 depends heavily on dominator. To use, you must first annotate a genbank file using domainate.py. For example:

`domainate.py -i input.gb -r reference.fasta --max_overlap 0.6 -o output.gb`

After annotating with a reference fasta file, you must generate two different enumerated reports. The first will be annotated by domain. For example:

`enum_report.py -i annotated_input.gb --by domain --domain_descriptions --start --end --score -o output_by_domain.tsv`

The second will be annotated by contig. For example:

`enum_report.py -i annotated_input.gb --by contig --length --topology -o output_by_contig.tsv`

After that, both reports can be fed into remt2. For example:

`python remt2.py -i by_domain.tsv  by_contig.tsv -o output`
*For inputs, by domain goes first and by contig goes second

remt2 is tailored for annotations using All_REBASE_Gold_Standards_Protein_20221213.fasta. To use a different reference file you will need to
get it in the same format as Gold Standard annotations. standardize.py (found in the helper programs folder) or remt2 could be tweaked to accomodate different 
fasta files. 

remt2 also contains a number of arguments. These include:
```
--coords will include the start and end points of the annotated domains listed
--stats prints a number of useful statistics about the run (thing like # of contigs, # of systems, etc.)
--excel makes the results excel friendly. The default puts the re/mt results on different lines and is better for reading HTML results
--filter_n will remove enzymes with recognition sequences that contain >3 Ns in a row. Can be useful as sequences with many Ns in a row can lead to some weird matches
--ingore_nonspec will remove enzyems with nonspecific recognition sequences (sequences <3 bp)
--min_distance takes an integer as an argument and will only output systems above the specified minimum distance
--min_score takes an integer as an argument and will only include enzymes above the specified score (scores come from domainate.py. Typically a min score of 40-60 is best)
--in_range takes an integer as an argument and will include the reference sequences of mts within the specified range of listed res and vice versa
```
