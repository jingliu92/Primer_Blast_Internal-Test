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
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq"  \
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
- alignment length ≥ 80% of primer length
- mismatches ≤ 4
```
2. For **Marker (gene)** presence, require:
```
- both F and R primers present
- same genome
- same contig (relax same-contig rule (draft assemblies) in the following try)
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

#Combine everything into one file**
awk 'NR==1 || $1!="Genome"' parsed/*.parsed.tsv > all_parsed.tsv
```
## Step 4: Collapse F/R → gene markers
**Pair primers + compute amplicon size**
` nano 2.all_parsed_advanced.py`

```
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
```
```
python all_parsed_advanced.py
```

