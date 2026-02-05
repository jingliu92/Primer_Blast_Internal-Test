#!/usr/bin/env python3
"""
collapse_marker_matrix.py

Convert marker_calls.tsv into a Genome × Marker
presence/absence matrix (1 = present, 0 = absent).
"""

import csv
from collections import defaultdict

INPUT = "marker_calls.tsv"
OUTPUT = "final_marker_matrix.tsv"

# Store presence information
presence = defaultdict(set)
genomes = set()
markers = set()

# Read marker calls
with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        genome = row["Genome"]
        marker = row["Marker"]

        genomes.add(genome)
        markers.add(marker)
        presence[genome].add(marker)

# Sort for stable output
genomes = sorted(genomes)
markers = sorted(markers)

# Write matrix
with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Genome"] + markers)

    for genome in genomes:
        row = [genome]
        for marker in markers:
            row.append(1 if marker in presence[genome] else 0)
        writer.writerow(row)

print(f"Wrote {OUTPUT}")
