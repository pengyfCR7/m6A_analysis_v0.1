#!/usr/bin/env bash
set -euo pipefail

input_bam=$1
output_bam=$2

mkdir -p "$(dirname "$output_bam")"

# Follow yuhaolab filterbam: retain alignments with -20 kb < TLEN < 20 kb.
# MAPQ30 and primary-alignment filtering were already applied to the input BAM.
samtools view -@ 4 -bh \
    -e 'tlen > -20000 && tlen < 20000' \
    -o "$output_bam" "$input_bam"

samtools index -@ 4 "$output_bam"
