 #!/usr/bin/env python3

from pathlib import Path
from itertools import product
from collections import Counter

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    ConfusionMatrixDisplay
)

project = Path(__file__).resolve().parent
data_folder = project / "data_clean"

species_labels = {
    "dmel": "Drosophila melanogaster",
    "dsim": "Drosophila simulans",
    "dsuz": "Drosophila suzukii"
}

# 3-mers
k = 3
all_kmers = [
    "".join(bases)
    for bases in product("ATGC", repeat=k)
]

# extract frequencies of 3-mers from a DNA sequence
def extract_kmers(sequence, k):
    sequence = sequence.upper()
    kmers = [
        sequence[i:i+k]
        for i in range(len(sequence) - k + 1)
    ]

    counts = Counter(kmers)

    total = len(kmers)
    frequencies = {
        mer: counts[mer] / total
        for mer in all_kmers
    }
    return frequencies

# process FASTA files
dataset = []
seen_sequences = {}
for fasta_file in sorted(data_folder.glob("*.fasta")):
    prefix = fasta_file.stem.split("_")[0]

    if prefix not in species_labels:
        continue
    species = species_labels[prefix]
    for record in SeqIO.parse(fasta_file, "fasta"):
        dna = str(record.seq).upper()

        if set(dna) - set("ATGC"):
            print("Skipping invalid sequence:", fasta_file.name)
            continue
        if len(dna) < k:
            continue

        if dna in seen_sequences:
            if seen_sequences[dna] != species:
                raise ValueError(
                    f"Conflicting species for sequence: {record.id}"
                )
            print("Skipping duplicate sequence:", record.id)
            continue
                
        seen_sequences[dna] = species

        features = extract_kmers(dna, k)
        features["accession"] = record.id
        features["filename"] = fasta_file.name
        features["species"] = species
        features["length_bp"] = len(dna)
        dataset.append(features)

df = pd.DataFrame(dataset)
if df.empty:
    raise ValueError("No valid sequences found in the dataset.")
print("\nSequences per species:")
print(df["species"].value_counts())

print("\nSequence lengths:")
print(df.groupby("species")["length_bp"].agg(
    ["min", "max", "mean"]
).round(1))

# Save the extracted k-mer features
df.to_csv(project / "kmer_features.csv", index=False)

# X contains the numerical k-mer features
X = df[all_kmers]

# y contains the correct species labels
y = df["species"]


# split training and dataset
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.25,
    random_state=42,
    stratify=y
)
print("\nTraining sequences:", len(X_train))
print("Testing sequences:", len(X_test))

# train the classifier
model = make_pipeline(
    StandardScaler(),
    LogisticRegression(
        max_iter=2000,
        class_weight="balanced",
    )
)

model.fit(X_train, y_train)

# Predict unseen sequences
predictions = model.predict(X_test)

print("\nMODEL PERFORMANCE")
print("-" * 40)
print("\nAccuracy:", accuracy_score(y_test, predictions))

print("\nClassification report:")
print(classification_report(
    y_test,
    predictions,
    zero_division=0
))


# baseline
baseline = DummyClassifier(strategy="most_frequent")
baseline.fit(X_train, y_train)

baseline_predictions = baseline.predict(X_test)

print("Baseline accuracy:", accuracy_score(y_test, baseline_predictions))


# confusion matrix

labels = list(species_labels.values())
ConfusionMatrixDisplay.from_predictions(
    y_test,
    predictions,
    labels=labels,
    cmap="Blues",
    values_format="d"
)

plt.title("Drosophila Species Classification (3-mers)")
plt.tight_layout()

plt.savefig(
    project / "confusion_matrix.png",
    dpi=300
)

plt.show()



from sklearn.model_selection import StratifiedKFold, cross_val_score

# Divide the dataset into five stratified folds
cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

# Evaluate the classifier across five folds
scores = cross_val_score(
    model,
    X,
    y,
    cv=cv,
    scoring="accuracy"
)

print("\nCROSS-VALIDATION RESULTS")
print("-" * 40)

for i, score in enumerate(scores):
    print(f"Fold {i+1}: {score:.2%}")

print(f"\nMean accuracy: {scores.mean():.2%}")
print(f"Standard deviation: {scores.std():.4f}")
