import os, sys
from Bio import SeqIO
INPUT_fasta = sys.argv[1]
OUTPUT_fasta = sys.argv[2]
SeqIO.convert(INPUT_fasta, 'fasta', OUTPUT_fasta, 'fasta')