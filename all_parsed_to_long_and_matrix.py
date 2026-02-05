#!/usr/bin/env python3
"""
all_parsed_to_long_and_matrix.py

From all_parsed.tsv:
1) Build marker_long.tsv (paired Fwd/Rev with amplicon size)
2) Build primer-level presence/absence matrix (all primers as columns)

No intermediate files required.
"""

import csv
from collections import defaultdict

INPUT = "all_parsed.tsv"

LONG_OUT = "marker_long.tsv"
MATRIX_OUT = "marker_matrix.tsv"

MIN_AMP = 100
MAX_AMP = 1500

# ----------------------------
# Read all_parsed.tsv
# ----------------------------
rows = []
genomes = set()
primers = set()

with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)
        genomes.add(row["Genome"])
        primers.add(row["Primer"])

genomes = sorted(genomes)
primers = sorted(primers)

# ----------------------------
# 1) Primer-level presence/absence matrix
# ----------------------------
primer_presence = defaultdict(set)

for r in rows:
    primer_presence[r["Genome"]].add(r["Primer"])

with open(MATRIX_OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Genome"] + primers)

    for genome in genomes:
        row = [genome]
        for primer in primers:
            row.append(1 if primer in primer_presence[genome] else 0)
        writer.writerow(row)

# ----------------------------
# 2) Marker-long table (Fwd + Rev pairing)
# ----------------------------
# Organize hits by genome + marker
hits = defaultdict(list)

for r in rows:
    primer = r["Primer"]

    if "Fwd" in primer:
        marker = primer.split("Fwd")[0]
        direction = "F"
    elif "Rev" in primer:
        marker = primer.split("Rev")[0]
        direction = "R"
    else:
        continue

    hits[(r["Genome"], marker)].append({
        "direction": direction,
        "primer": primer,
        "contig": r["Contig"],
        "strand": r["Strand"],
        "start": int(r["Start"]),
        "primer_len": r["Primer_Length"],
        "align_len": r["Align_Length"],
        "pid": r["Percent_Identity"],
        "qseq": r["Query_Align_Seq"],
        "sseq": r["Subject_Align_Seq"]
    })

long_rows = []

for (genome, marker), records in hits.items():
    fwds = [x for x in records if x["direction"] == "F"]
    revs = [x for x in records if x["direction"] == "R"]

    for f in fwds:
        for r in revs:
            # opposite strand required
            if f["strand"] == r["strand"]:
                continue

            amp = abs(f["start"] - r["start"])
            if not (MIN_AMP <= amp <= MAX_AMP):
                continue

            long_rows.append([
                genome,
                marker,
                f["primer"],
                r["primer"],
                f["contig"],
                r["contig"],
                amp,
                f["primer_len"],
                f["align_len"],
                f["pid"],
                f["qseq"],
                f["sseq"],
                r["primer_len"],
                r["align_len"],
                r["pid"],
                r["qseq"],
                r["sseq"]
            ])

# Write marker_long.tsv
with open(LONG_OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Genome",
        "Marker",
        "Forward_Primer",
        "Reverse_Primer",
        "Forward_Contig",
        "Reverse_Contig",
        "Amplicon_Size",
        "F_Primer_Length",
        "F_Align_Length",
        "F_Percent_Identity",
        "F_Query_Align_Seq",
        "F_Subject_Align_Seq",
        "R_Primer_Length",
        "R_Align_Length",
        "R_Percent_Identity",
        "R_Query_Align_Seq",
        "R_Subject_Align_Seq"
    ])
    writer.writerows(long_rows)

print("Generated:")
print(f"  {LONG_OUT}")
print(f"  {MATRIX_OUT}")
