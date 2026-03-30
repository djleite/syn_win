# syn_win
Quick and dirty chatGPT synteny tool. Have reference gene, plus up/down stream genes, and blast to other species. Terrible plot, and needs checking for sanity.

## installation with install_env.sh
you can conda install python and diamond
pip install dependencies

## data
I used GCF_000001215.4 and GCF_030788295.1 Drosophila genomes to test. Only need protein and GFF files.

## Run

'python pipeline.py'

You need to to change the script at the bottom to iput the reference and test species and the gene name of the reference. IMPORTANT just use the first name in the protein fasta ID.

## Output
- HTML interactive
- PDF
- SVG

