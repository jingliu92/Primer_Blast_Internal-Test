# Primer_Blast_Internal-Test

We have 1,353 isolates in our collection. 

## Step 1: Loop through genomes and BLAST primers
```
mkdir -p dbs results

for g in /home/jing/E.coli_test/ecoli_all/*.contigs.fasta
do
  genome=$(basename "$g" .contigs.fasta)

  echo "Processing $genome"

  makeblastdb \
    -in "$g" \
    -dbtype nucl \
    -out dbs/$genome \
    > /dev/null

  blastn \
    -query primer_markers.fasta \
    -db dbs/$genome \
    -task blastn-short \
    -word_size 7 \
    -evalue 1000 \
    -dust no \
    -soft_masking false \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -out results/${genome}.tsv
done
```
#### Example Results {genome}.tsv 
<img width="1012" height="221" alt="image" src="https://github.com/user-attachments/assets/57696906-8de1-4350-b075-b115fa3c569f" />
We notice that there are so many hits for each of the primer, this is because our primer is less than 30bp. That’s short.

**BLAST behavior for short sequences:**
- finds partial matches
- finds internal fragments
- finds low-complexity k-mer hits
- finds matches on both strands
- repeats this across every contig

So these lines are NOT full primer binding sites and biologically meaningless for PCR. We must ONLY keep near-full-length matches.

## Step 2: Define biological rules
1. For **primer presence**, require:
```
- alignment length ≥ 90% of primer length
- mismatches ≤ 2
```
2. For **Marker (gene)** presence, require:
```
- both F and R primers present
- same genome
- same contig
- opposite orientation
- amplicon size within range
```
## Step 3: Parse BLAST Result → primer hits (Python)
`parse_primer_hits.py`
```
#!/usr/bin/env python3
"""
parse_primer_hits.py

Parse BLAST tabular output for primer searches and retain
only biologically meaningful primer binding sites.

Expected BLAST outfmt:
6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore

Filters:
- alignment length >= 90% of primer length
- mismatches <= 2   (approximated as qlen - length)

Outputs:
Genome, Primer, Contig, Strand,
Start, Stop,
Primer_Length, Align_Length,
Percent_Identity
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

        # -------- STRICT FILTERS --------
        coverage = length / qlen
        mismatches = qlen - length

        if coverage < 0.90:
            continue
        if mismatches > 10:
            continue
        # --------------------------------

        strand = "+" if sstart < send else "-"
        start  = min(sstart, send)
        stop   = max(sstart, send)

        hits[qseqid].append(
            (sseqid, strand, start, stop, qlen, length, pident)
        )

# Print header (report-friendly)
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
    sep="\t"
)

for primer, records in hits.items():
    for contig, strand, start, stop, qlen, length, pident in records:
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
            sep="\t"
        )
```
**Batch Run**
```
mkdir -p parsed

for f in results/*.tsv
do
  genome=$(basename "$f" .tsv)

  echo "Parsing $genome"

  python 1.parse_primer_hits.py "$f" "$genome" > parsed/${genome}.parsed.tsv
done
```
**Combine everything into one file**
```
awk 'NR==1 || $1!="Genome"' parsed/*.parsed.tsv > all_parsed.tsv
```
```
