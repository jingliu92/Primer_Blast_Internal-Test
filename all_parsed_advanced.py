#!/usr/bin/env python3
"""
all_parsed_advanced.py

Advanced primer pairing from all_parsed.tsv with:
- relaxed same-contig rule (draft assemblies)
- amplicon size metadata
- primer length, alignment length, % identity
- alignment DNA sequences (query + subject)
- long-format + wide-format output
- multiple primer variants per gene
- QC warnings for multiple loci
"""

import csv
from collections import defaultdict

INPUT = "all_parsed.tsv"

# Outputs
LONG_OUT = "marker_long.tsv"
MATRIX_OUT = "marker_matrix.tsv"
QC_OUT = "marker_qc.tsv"

# Amplicon constraints
MIN_AMP = 100
MAX_AMP = 1500

# Store primer hits
hits = defaultdict(list)
genomes = set()
markers = set()

# ----------------------------
# Read all_parsed.tsv
# ----------------------------
with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        genome = row["Genome"]
        primer = row["Primer"]
        contig = row["Contig"]
        strand = row["Strand"]
        start = int(row["Start"])

        primer_len = int(row["Primer_Length"])
        align_len = int(row["Align_Length"])
        pid = float(row["Percent_Identity"])

        qseq = row["Query_Align_Seq"]
        sseq = row["Subject_Align_Seq"]

        genomes.add(genome)

        if "Fwd" in primer:
            marker = primer.split("Fwd")[0]
            direction = "F"
        elif "Rev" in primer:
            marker = primer.split("Rev")[0]
            direction = "R"
        else:
            continue

        markers.add(marker)

        hits[(genome, marker)].append({
            "direction": direction,
            "primer": primer,
            "contig": contig,
            "strand": strand,
            "start": start,
            "primer_len": primer_len,
            "align_len": align_len,
            "pid": pid,
            "qseq": qseq,
            "sseq": sseq
        })

# ----------------------------
# Pair primers & collect data
# ----------------------------
long_rows = []
presence = defaultdict(set)
qc_warnings = []

for (genome, marker), records in hits.items():
    fwds = [r for r in records if r["direction"] == "F"]
    revs = [r for r in records if r["direction"] == "R"]

    loci_found = []

    for f in fwds:
        for r in revs:
            # must be opposite strand
            if f["strand"] == r["strand"]:
                continue

            amp = abs(f["start"] - r["start"])
            if not (MIN_AMP <= amp <= MAX_AMP):
                continue

            loci_found.append((f, r, amp))

    if not loci_found:
        continue

    # Marker is present
    presence[genome].add(marker)

    # QC: multiple loci
    if len(loci_found) > 1:
        qc_warnings.append([
            genome, marker, "MULTIPLE_LOCI", len(loci_found)
        ])

    # Long-format rows
    for f, r, amp in loci_found:
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

# ----------------------------
# Write long-format table
# ----------------------------
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

# ----------------------------
# Write wide presence/absence matrix
# ----------------------------
genomes = sorted(genomes)
markers = sorted(markers)

with open(MATRIX_OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Genome"] + markers)
    for genome in genomes:
        row = [genome]
        for marker in markers:
            row.append(1 if marker in presence[genome] else 0)
        writer.writerow(row)

# ----------------------------
# Write QC warnings
# ----------------------------
with open(QC_OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Genome", "Marker", "Issue", "Count"])
    writer.writerows(qc_warnings)

print("Generated:")
print(f"  {LONG_OUT}")
print(f"  {MATRIX_OUT}")
print(f"  {QC_OUT}")
