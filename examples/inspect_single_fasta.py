#!/usr/bin/env python3


from Bio import SeqIO
from collections import Counter

# Read downloaded DNA sequences from a FASTA file
record = SeqIO.read("data/dmel_01.fasta", "fasta")

# convert sequence to uppercase string
dna = str(record.seq).upper()

# count each nucleotide
counts = Counter(dna)

# calculate GC content
gc_content = (counts['G'] + counts['C']) / len(dna) * 100

print("Accession:", record.id)
print("Description:", record.description)
print("Sequence Length:", len(dna), "bp")
print("Nucleotide Counts:", dict(counts))
print("GC Content:", f"{gc_content:.2f} %")