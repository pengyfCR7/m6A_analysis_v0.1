#!/usr/bin/env bash
set -euo pipefail

project=/data/home/pengyf/projects/m6A/Defining_context_dependent_m6A_PRJNA1028096
input_dir="$project/data/fastp_STAR_TAIR10_Araport11_current/star_mapq30_primary"
output_dir="$project/data/yuhaolab_MAPQ30_PRIMARY_SECOND_STRAND"

source /etc/profile.d/lmod.sh
module purge
module load samtools/1.24

mkdir -p "$output_dir/20k" "$output_dir/stranded"
for input_bam in "$input_dir"/*.bam; do
    name=$(basename "$input_bam" .bam)
    bam20k="$output_dir/20k/${name}.20k.bam"
    echo "[$(date '+%F %T')] $name"
    bash "$project/scripts02/filter_20k_bam.sh" "$input_bam" "$bam20k"
    bash "$project/scripts02/split_bam_second_strand.sh" "$bam20k" "$output_dir/stranded"
done
