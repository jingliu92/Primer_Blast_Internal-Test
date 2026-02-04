# Primer_Blast_Internal-Test

We have 1,353 isolates in our collection. 

Step 1: Loop through genomes and BLAST primers
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
So these lines are NOT full primer binding sites and biologically meaningless for PCR.
STEP 2. Define biological rules


Primer hit is valid if:

length / qlen ≥ 0.9

mismatch ≤ 2

Marker (gene) is present if:

both F and R primers present

same genome

same contig

opposite orientation

amplicon size within range

Example amplicon size:
Step 2:Parse BLAST Result → primer hits (Python)
