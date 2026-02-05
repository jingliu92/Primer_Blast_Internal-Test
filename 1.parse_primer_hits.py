#!/usr/bin/env python3
"""
parse_primer_hits.py

Parse BLAST tabular output for primer searches and retain
biologically meaningful primer binding sites, including
alignment sequences.

Expected BLAST outfmt:
6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq

Filters:
- alignment length >= 90% of primer length
- mismatches <= 3   (approximated as qlen - length)
- percent identity >= 30%

Outputs:
Genome
Primer
Contig
Strand
Start
Stop
Primer_Length
Align_Length
Percent_Identity
Query_Align_Seq
Subject_Align_Seq
"""

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    sys.stderr.write(
        "Usage: parse_primer_hits.py <blast.tsv> <genome_id>\n"
    )
    sys.exit(1)

blast_file = sys.argv[1]
genome_id = sys.argv[2]

hits = defaultdict(list)

with open(blast_file) as fh:
    for line in fh:
        if not line.strip():
            continue

        cols = line.rstrip().split("\t")

        qseqid = cols[0]
        sseqid = cols[1]
        pident = float(cols[2])
        length = int(cols[3])
        qlen   = int(cols[4])
        sstart = int(cols[7])
        send   = int(cols[8])
        qseq   = cols[11]
        sseq   = cols[12]

        # -------- FILTERS --------
        coverage = length / qlen
        mismatches = qlen - length

        if coverage < 0.90:
            continue
        if mismatches > 3:
            continue
        if pident < 30.0:
            continue
        # -------------------------

        strand = "+" if sstart < send else "-"
        start  = min(sstart, send)
        stop   = max(sstart, send)

        hits[qseqid].append(
            (
                sseqid,
                strand,
                start,
                stop,
                qlen,
                length,
                pident,
                qseq,
                sseq
            )
        )

# Header
print(
    "Genome",
    "Primer",
    "Contig",
    "Strand",
    "Start",
    "Stop",
    "Primer_Length",
    "Align_Length",
    "Percent_Identity",
    "Query_Align_Seq",
    "Subject_Align_Seq",
    sep="\t"
)

# Output records
for primer, records in hits.items():
    for (
        contig,
        strand,
        start,
        stop,
        qlen,
        length,
        pident,
        qseq,
        sseq
    ) in records:
        print(
            genome_id,
            primer,
            contig,
            strand,
            start,
            stop,
            qlen,
            length,
            f"{pident:.2f}",
            qseq,
            sseq,
            sep="\t"
        )
