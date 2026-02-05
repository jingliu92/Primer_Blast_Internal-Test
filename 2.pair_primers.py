#!/usr/bin/env python3
"""
pair_primers.py

From all_parsed.tsv:
- Pair Fwd + Rev primers
- Require same genome
- Require same contig
- Require opposite strand
- Compute amplicon size
- Call marker presence
"""

from collections import defaultdict

INPUT = "all_parsed.tsv"
MIN_AMP = 100     # adjust if needed
MAX_AMP = 1500    # adjust if needed

hits = defaultdict(list)

# Read input
with open(INPUT) as f:
    header = f.readline()
    for line in f:
        genome, primer, contig, strand, start, stop, qlen, alen, pid = line.strip().split("\t")
        start, stop = int(start), int(stop)

        # derive marker + direction
        if "Fwd" in primer:
            marker = primer.split("Fwd")[0]
            direction = "F"
        elif "Rev" in primer:
            marker = primer.split("Rev")[0]
            direction = "R"
        else:
            continue

        hits[(genome, marker)].append(
            (direction, contig, strand, start, stop)
        )

# Output
print("Genome", "Marker", "Contig", "Amplicon_Size", sep="\t")

for (genome, marker), records in hits.items():
    fwds = [r for r in records if r[0] == "F"]
    revs = [r for r in records if r[0] == "R"]

    for f in fwds:
        for r in revs:
            _, contig_f, strand_f, start_f, stop_f = f
            _, contig_r, strand_r, start_r, stop_r = r

            # same contig
            if contig_f != contig_r:
                continue

            # opposite strand
            if strand_f == strand_r:
                continue

            # amplicon size
            amp = abs(start_f - start_r)

            if MIN_AMP <= amp <= MAX_AMP:
                print(genome, marker, contig_f, amp, sep="\t")
                break
