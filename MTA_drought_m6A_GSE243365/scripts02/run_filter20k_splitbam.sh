#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

project=/data/home/pengyf/projects/m6A/MTA_drought_m6A_GSE243365
input_dir="$project/data/fastp_STAR_TAIR10_Araport11_current/star_mapq30_primary"
output_dir="$project/data/yuhaolab_MAPQ30_PRIMARY"
scripts="$project/scripts02"

source /etc/profile.d/lmod.sh
module purge
module load samtools/1.24

inputs=("$input_dir"/{WT,mta}_{Unstressed,Drought}_{m6A_IP,Input}_R*.star.unique.bam)
if (( ${#inputs[@]} != 24 )); then
    echo "ERROR: expected 24 m6A IP/Input BAMs; found ${#inputs[@]}" >&2
    exit 1
fi

mkdir -p "$output_dir/20k" "$output_dir/stranded"
for input_bam in "${inputs[@]}"; do
    name=$(basename "$input_bam" .bam)
    bam20k="$output_dir/20k/${name}.20k.bam"
    echo "[$(date '+%F %T')] $name"
    bash "$scripts/filter_20k_bam.sh" "$input_bam" "$bam20k"
    bash "$scripts/split_bam_reverse_stranded.sh" "$bam20k" "$output_dir/stranded"
done
