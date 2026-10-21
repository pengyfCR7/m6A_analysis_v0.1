#!/usr/bin/env bash
set -euo pipefail

input_bam=$1
output_dir=$2

mkdir -p "$output_dir"
name=$(basename "$input_bam" .bam)

# RSeQC: m6A IP and Input are reverse/first-strand; Read 1 is antisense to RNA.
samtools view -@ 4 -bh \
    -e '((flag & 83) == 83) || ((flag & 163) == 163)' \
    -o "$output_dir/${name}.fwd.bam" "$input_bam"

samtools view -@ 4 -bh \
    -e '((flag & 99) == 99) || ((flag & 147) == 147)' \
    -o "$output_dir/${name}.rev.bam" "$input_bam"

samtools index -@ 4 "$output_dir/${name}.fwd.bam"
samtools index -@ 4 "$output_dir/${name}.rev.bam"
