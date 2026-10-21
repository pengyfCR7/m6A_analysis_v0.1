#!/usr/bin/env bash
set -euo pipefail

project=/data/home/pengyf/projects/m6A/Defining_context_dependent_m6A_PRJNA1028096
base="$project/data/yuhaolab_MAPQ30_PRIMARY_SECOND_STRAND"
sample_table="$project/config/sample_exomepeak2_temperature_and_light_quality_mapq30_primary.txt"

source /etc/profile.d/lmod.sh
module purge
module load MACS3/3.0.4

mkdir -p "$base/macs3_stranded"
awk -F '\t' 'NR > 1 && $2 ~ /_IP_/ {print $2}' "$sample_table" | while read -r sample; do
    control=${sample/_IP_/_input_}
    for direction in fwd rev; do
        name="$sample.$direction"
        out="$base/macs3_stranded/$name"
        mkdir -p "$out"
        echo "[$(date '+%F %T')] $name"
        macs3 callpeak \
            -t "$base/stranded/$sample.star.unique.20k.$direction.bam" \
            -c "$base/stranded/$control.star.unique.20k.$direction.bam" \
            -f BAMPE -g 120000000 -n "$name" -B \
            --nomodel --call-summits --min-length 100 \
            --keep-dup all -q 0.05 --outdir "$out"
    done
done
