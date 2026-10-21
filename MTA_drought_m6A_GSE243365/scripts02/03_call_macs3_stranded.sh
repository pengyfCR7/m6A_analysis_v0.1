#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

source /etc/profile.d/lmod.sh
module purge
module load MACS3/3.0.4

project=/data/home/pengyf/projects/m6A/MTA_drought_m6A_GSE243365
bam_dir="$project/data/yuhaolab_MAPQ30_PRIMARY/stranded"
out_dir="$project/data/yuhaolab_MAPQ30_PRIMARY/macs3_stranded"

ips=("$bam_dir"/{WT,mta}_{Unstressed,Drought}_m6A_IP_R[123].star.unique.20k.{fwd,rev}.bam)
if (( ${#ips[@]} != 24 )); then
    echo "ERROR: expected 24 stranded IP BAMs; found ${#ips[@]}" >&2
    exit 1
fi

mkdir -p "$out_dir"
macs3 --version > "$out_dir/macs3_version.txt" 2>&1
printf 'sample\tstrand\tpeak_rows\tsummit_rows\tip_bam\tinput_bam\n' > "$out_dir/peak_counts.tsv"

for ip in "${ips[@]}"; do
    input=${ip/_m6A_IP_/_Input_}
    base=$(basename "$ip")
    stem=${base%.star.unique.20k.*.bam}
    strand=${base%.bam}; strand=${strand##*.}
    sample=${stem/_m6A_IP_/_}
    name="$sample.$strand"
    sample_out="$out_dir/$name"

    [[ -s "$input" ]] || { echo "ERROR: missing Input: $input" >&2; exit 1; }
    mkdir -p "$sample_out"
    echo "[CALL] $name"
    macs3 callpeak \
        -t "$ip" -c "$input" \
        -f BAMPE -g 120000000 \
        --nomodel --call-summits --min-length 100 --keep-dup all \
        -q 0.05 -n "$name" --outdir "$sample_out" \
        > "$sample_out/macs3.log" 2>&1

    peaks="$sample_out/${name}_peaks.narrowPeak"
    summits="$sample_out/${name}_summits.bed"
    [[ -f "$peaks" ]] || { echo "ERROR: missing $peaks" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$strand" "$(wc -l < "$peaks")" "$(wc -l < "$summits")" "$ip" "$input" \
        >> "$out_dir/peak_counts.tsv"
done

echo "Summary: $out_dir/peak_counts.tsv"
