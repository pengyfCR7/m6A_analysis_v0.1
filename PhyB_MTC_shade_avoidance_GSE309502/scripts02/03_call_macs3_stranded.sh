#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

source /etc/profile.d/lmod.sh
module purge
module load MACS3/3.0.4
# GSE309502: stranded BAMs already represent RNA transcript strand.
ROOT="${ROOT:-/data/home/pengyf/projects/m6A/PhyB_MTC_shade_avoidance_GSE309502}"
BAM_DIR="$ROOT/data/yuhaolab_MAPQ30_PRIMARY/stranded"
OUT_DIR="$ROOT/data/yuhaolab_MAPQ30_PRIMARY/macs3_stranded"
DRY_RUN="${DRY_RUN:-0}"

# Each IP library has a matching Input with the same condition, replicate and strand.
ips=("$BAM_DIR"/*_IP_R[12].star.unique.20k.fwd.bam
     "$BAM_DIR"/*_IP_R[12].star.unique.20k.rev.bam)

if (( ${#ips[@]} != 12 )); then
    echo "ERROR: Expected 12 IP BAMs; found ${#ips[@]} in $BAM_DIR" >&2
    exit 1
fi

for ip in "${ips[@]}"; do
    input="${ip/_IP_/_Input_}"
    if [[ ! -s "$input" ]]; then
        echo "ERROR: Matching Input BAM missing: $input" >&2
        exit 1
    fi
done

if [[ "$DRY_RUN" != 1 ]]; then
    command -v macs3 >/dev/null || { echo 'ERROR: macs3 not found' >&2; exit 1; }
    mkdir -p "$OUT_DIR"
    macs3 --version > "$OUT_DIR/macs3_version.txt" 2>&1
    printf 'sample\tstrand\tpeak_rows\tsummit_rows\tip_bam\tinput_bam\n' > "$OUT_DIR/peak_counts.tsv"
fi

for ip in "${ips[@]}"; do
    input="${ip/_IP_/_Input_}"
    base="$(basename "$ip")"
    stem="${base%.star.unique.20k.*.bam}"
    strand="${base%.bam}"; strand="${strand##*.}"
    sample="${stem/_IP_/_}"
    name="${sample}.${strand}"
    sample_out="$OUT_DIR/$name"
    peaks="$sample_out/${name}_peaks.narrowPeak"
    summits="$sample_out/${name}_summits.bed"

    if [[ "$DRY_RUN" == 1 ]]; then
        printf '[DRY RUN] %-28s %s\n          Input: %s\n' "$name" "$ip" "$input"
        continue
    fi

    mkdir -p "$sample_out"
    if [[ -f "$peaks" ]]; then
        echo "[SKIP] $name (narrowPeak already exists)"
    else
        echo "[CALL] $name"
        macs3 callpeak \
            -t "$ip" -c "$input" \
            -f BAMPE -g 120000000 \
            --nomodel --call-summits --min-length 100 --keep-dup all \
            -q 0.05 \
            -n "$name" --outdir "$sample_out" \
            > "$sample_out/macs3.log" 2>&1
    fi

    [[ -f "$peaks" ]] || { echo "ERROR: Missing MACS3 narrowPeak: $peaks" >&2; exit 1; }
    peak_n=$(wc -l < "$peaks")
    summit_n=0
    if [[ -f "$summits" ]]; then summit_n=$(wc -l < "$summits"); fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$strand" "$peak_n" "$summit_n" "$ip" "$input" \
        >> "$OUT_DIR/peak_counts.tsv"
    printf '[DONE] %-28s peaks=%s  summits=%s\n' "$name" "$peak_n" "$summit_n"
done

if [[ "$DRY_RUN" != 1 ]]; then
    echo "Summary: $OUT_DIR/peak_counts.tsv"
fi
