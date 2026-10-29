#!/bin/sh
# Validation, step 2 (docs/VALIDATION.md): BLAST of the queries against the
# reference, with wide limits, as an independent source of candidate pairs.
#   sh tools/validation/run_blast.sh REFERENCE.fasta WORKDIR
# Needs makeblastdb and blastn (BLAST+ 2.17.0 was used) on the PATH.
set -e
ref="$1"; work="$2"
mkdir -p "$work/blastdb"
cp "$ref" "$work/blastdb/ref.fasta"
makeblastdb -in "$work/blastdb/ref.fasta" -dbtype nucl -out "$work/blastdb/ref" > /dev/null
blastn -task blastn -query "$work/queries.fasta" -db "$work/blastdb/ref" -word_size 11   -perc_identity 80 -qcov_hsp_perc 70 -max_target_seqs 20000 -max_hsps 1 -evalue 1e-10   -dust no -soft_masking false -strand both -num_threads 4   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand qlen"   > "$work/blast_hits.tsv"
