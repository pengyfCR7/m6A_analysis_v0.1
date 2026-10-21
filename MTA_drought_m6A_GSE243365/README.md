# MTA drought m⁶A：GSE243365

该项目按reverse/first-strand处理，包括WT/mta、Unstressed/Drought各3个重复。

```bash
bash scripts02/run_filter_split_macs3.sh
THREADS=8 /opt/conda-envs/exomePeak2/1.14.3/bin/Rscript scripts02/run_m6A_union_poisson_compact.R
```

实际结果：`results/stranded_m6A_union_poisson_compact/union_poisson_all_peaks.tsv`。
