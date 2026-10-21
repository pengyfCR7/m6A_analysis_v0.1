#!/usr/bin/env bash
set -euo pipefail

project=/data/home/pengyf/projects/m6A/MTA_drought_m6A_GSE243365
bash "$project/scripts02/run_filter20k_splitbam.sh"
bash "$project/scripts02/03_call_macs3_stranded.sh"
