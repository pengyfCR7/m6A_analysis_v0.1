#!/usr/bin/env bash
set -euo pipefail

project=/data/home/pengyf/projects/m6A/Defining_context_dependent_m6A_PRJNA1028096
cd "$project"

bash scripts02/01_filter20k_splitbam_second_strand.sh
bash scripts02/02_call_macs3_stranded.sh

THREADS=8 /opt/conda-envs/exomePeak2/1.14.3/bin/Rscript \
    scripts02/run_m6A_union_poisson_compact.R
