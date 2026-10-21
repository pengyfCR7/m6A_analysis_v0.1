# PhyB/MTC shade：GSE309502

该项目经RSeQC确认按reverse/first-strand处理，`split_bam_reverse_stranded.sh`中的fwd/rev规则与second-strand项目相反。

```bash
bash scripts02/run_filter20k_splitbam.sh
bash scripts02/03_call_macs3_stranded.sh
THREADS=8 /opt/conda-envs/exomePeak2/1.14.3/bin/Rscript scripts02/run_m6A_union_poisson_compact.R
```

实际结果：`results/stranded_m6A_union_poisson_compact/union_poisson_all_peaks.tsv`。
