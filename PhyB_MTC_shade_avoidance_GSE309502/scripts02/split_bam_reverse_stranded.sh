#!/usr/bin/env bash
set -euo pipefail

input_bam=$1
output_dir=$2

mkdir -p "$output_dir"
name=$(basename "$input_bam" .bam)
fwd_bam="$output_dir/${name}.fwd.bam"
rev_bam="$output_dir/${name}.rev.bam"

# GSE309502 is reverse-stranded: Read 1 is antisense to the RNA.
# Therefore fwd/rev are the reverse of yuhaolab's Read-1-sense splitBam labels.
samtools view -@ 4 -bh \
    -e '((flag & 83) == 83) || ((flag & 163) == 163)' \
    -o "$fwd_bam" "$input_bam"

samtools view -@ 4 -bh \
    -e '((flag & 99) == 99) || ((flag & 147) == 147)' \
    -o "$rev_bam" "$input_bam"

samtools index -@ 4 "$fwd_bam"
samtools index -@ 4 "$rev_bam"
