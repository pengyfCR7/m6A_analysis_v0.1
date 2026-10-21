#!/usr/bin/env bash
set -euo pipefail

input_bam=$1
output_bam=$2

mkdir -p "$(dirname "$output_bam")"
samtools view -@ 4 -bh \
    -e 'tlen > -20000 && tlen < 20000' \
    -o "$output_bam" "$input_bam"
samtools index -@ 4 "$output_bam"
