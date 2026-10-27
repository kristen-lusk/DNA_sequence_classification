 #!/usr/bin/env python3
 

from pathlib import Path
from Bio import SeqIO
from collections import Counter
import csv

# Locate the data folder
data_folder = Path(__file__).resolve().parent / "data"

# Match filenames to species
species_labels = {
    "dmel": "Drosophila melanogaster",
    "dsim": "Drosophila simulans",
    "dsuz": "Drosophila suzukii"
}

results = []
seen_sequences = {}

# Loop through every FASTA file
for fasta_file in sorted(data_folder.glob("*.fasta")):

    prefix = fasta_file.stem.split("_")[0]

    if prefix not in species_labels:
        print("Skipping unrecognized file:", fasta_file.name)
        continue

    species = species_labels[prefix]

    # Read every sequence within the file
    for record in SeqIO.parse(fasta_file, "fasta"):

        dna = str(record.seq).upper()
        length = len(dna)

        if length == 0:
            print("Empty sequence:", fasta_file.name)
            continue

        # Calculate nucleotide frequencies
        counts = Counter(dna)

        gc_content = (
            counts["G"] + counts["C"]
        ) / length * 100

        # Check for unexpected characters
        invalid_bases = set(dna) - set("ATGC")

        # Check for duplicate DNA sequences
        duplicate_of = seen_sequences.get(dna, "")

        if not duplicate_of:
            seen_sequences[dna] = record.id

        # Check whether the description matches the label
        species_match = (
            species.lower() in record.description.lower()
        )

        results.append({
            "filename": fasta_file.name,
            "accession": record.id,
            "species": species,
            "length_bp": length,
            "gc_percent": round(gc_content, 2),
            "species_match": species_match,
            "invalid_bases": "".join(sorted(invalid_bases)),
            "duplicate_of": duplicate_of
        })

# Summarize the dataset
print("\nDNA SEQUENCE DATASET SUMMARY")
print("-" * 40)

for species in species_labels.values():

    subset = [
        r for r in results
        if r["species"] == species
    ]

    if not subset:
        continue

    lengths = [r["length_bp"] for r in subset]

    print(f"\nSpecies: {species}")
    print(f"Sequences: {len(subset)}")
    print(f"Minimum length: {min(lengths)} bp")
    print(f"Maximum length: {max(lengths)} bp")

print("\nTotal sequences:", len(results))

print("Exact duplicates:", sum(
    bool(r["duplicate_of"]) for r in results
))

print("Species label mismatches:", sum(
    not r["species_match"] for r in results
))

print("Sequences with invalid bases:", sum(
    bool(r["invalid_bases"]) for r in results
))

# Save results to a CSV file
output_file = data_folder.parent / "sequence_quality_report.csv"

with open(output_file, "w", newline="") as file:

    writer = csv.DictWriter(
        file,
        fieldnames=results[0].keys()
    )

    writer.writeheader()
    writer.writerows(results)

print("\nReport saved:", output_file)




# Find the project directory
project = Path(__file__).resolve().parent

# Input and output paths
report_file = project / "sequence_quality_report.csv"
raw_folder = project / "data"
clean_folder = project / "data_clean"

clean_folder.mkdir(exist_ok=True)

# Conflicting species identification:
# dmel_13 and dsim_02 contain identical DNA
# despite having different species labels.
flagged_files = {"dmel_13.fasta", "dsim_02.fasta"}

kept = Counter()
excluded = Counter()

with open(report_file, newline="") as file:
    reader = csv.DictReader(file)

    for row in reader:
        filename = row["filename"]
        species = row["species"]

        # Remove exact duplicates
        if row["duplicate_of"]:
            excluded["duplicate"] += 1
            continue

        # Exclude conflicting species identification
        if filename in flagged_files:
            excluded["species_conflict"] += 1
            continue

        # Exclude sequences with ambiguous nucleotides
        if row["invalid_bases"]:
            excluded["ambiguous_bases"] += 1
            continue

        # Verify the FASTA exists
        input_file = raw_folder / filename

        if not input_file.exists():
            print("File not found:", filename)
            continue

        # Read and save the accepted sequence
        record = SeqIO.read(input_file, "fasta")

        output_file = clean_folder / filename
        SeqIO.write(record, output_file, "fasta")

        kept[species] += 1

print("\nCLEANED DATASET")
print("-" * 40)

for species, count in kept.items():
    print(f"{species}: {count}")

print("\nTotal retained:", sum(kept.values()))
print("Exclusion reasons:", dict(excluded))
print("Saved to:", clean_folder)
