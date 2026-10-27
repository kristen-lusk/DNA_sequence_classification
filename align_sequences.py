 #!/usr/bin/env python3
 

from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import csv

# Locate cleaned sequences
project = Path(__file__).resolve().parent
data_folder = project / "data_clean"

# Load reference COI sequence
reference = SeqIO.read(
    data_folder / "dmel_01.fasta", "fasta"
)
ref_dna = reference.seq.upper()

# Configure local sequence alignment
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

results = []

# Compare every sequence with reference
for fasta_file in sorted(data_folder.glob("*.fasta")):

    for record in SeqIO.parse(fasta_file, "fasta"):
        dna = record.seq.upper()

        # Check both DNA orientations
        forward = aligner.align(ref_dna, dna)[0]
        reverse = aligner.align(
            ref_dna, dna.reverse_complement()
        )[0]

        best = max(
            [forward, reverse],
            key=lambda alignment: alignment.score
        )

        orientation = (
            "forward" if forward.score >= reverse.score
            else "reverse"
        )

        query = (
            dna if orientation == "forward"
            else dna.reverse_complement()
        )

        # Identify aligned nucleotide positions
        pairs = [
            (i, j)
            for i, j in best.indices.T
            if i >= 0 and j >= 0
        ]

        aligned_bp = len(pairs)

        matches = sum(
            ref_dna[i] == query[j]
            for i, j in pairs
        )

        identity = 100 * matches / aligned_bp
        coverage = 100 * aligned_bp / len(ref_dna)

        ref_start = min(i for i, j in pairs) + 1
        ref_end = max(i for i, j in pairs) + 1

        results.append({
            "filename": fasta_file.name,
            "accession": record.id,
            "aligned_bp": aligned_bp,
            "identity_percent": round(identity, 2),
            "reference_coverage_percent": round(
                coverage, 2
            ),
            "reference_start": ref_start,
            "reference_end": ref_end,
            "orientation": orientation
        })

# Save alignment report
output = project / "alignment_quality_report.csv"

with open(output, "w", newline="") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=results[0].keys()
    )
    writer.writeheader()
    writer.writerows(results)

# 6. Show sequences needing closer inspection
print("Sequences analyzed:", len(results))

print("\nSequences to review:")
for row in results:
    if (row["reference_coverage_percent"] < 75
            or row["identity_percent"] < 75):
        print(
            row["filename"],
            "| Identity:", row["identity_percent"],
            "| Coverage:",
            row["reference_coverage_percent"]
        )

print("\nReport saved to:", output)
