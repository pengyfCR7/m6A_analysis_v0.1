# m6A_analysis_v0.1

目前正在使用的、以链特异性 MACS3 union peaks 为基础的 m⁶A-seq 探索性分析脚本。仓库保留三个已经实际跑通的项目版本，不把不同文库方向和不同实验设计强行包装成一个复杂的通用流程。

## 已运行项目

- `PhyB_MTC_shade_avoidance_GSE309502`：shade、phyB；reverse/first-strand。
- `MTA_drought_m6A_GSE243365`：WT/mta drought；reverse/first-strand。
- `Defining_context_dependent_m6A_PRJNA1028096`：温度和光质；second-strand。

每个项目的基本步骤相同：

```text
MAPQ30 primary paired-end BAM
→ 保留 -20 kb < TLEN < 20 kb
→ 按实际RNA链拆成 fwd/rev BAM
→ 每个IP/Input样品、每条链分别运行MACS3
→ 同链MACS3 peaks取union
→ featureCounts统计union peak和代表转录本exon fragments
→ 以exon counts估计IP/Input各自的size factors
→ 配对Poisson Wald GLM计算m⁶A变化
```

依赖：`samtools 1.24`、`MACS3 3.0.4`、R及`GenomicRanges`、`Rsubread`、`DESeq2`。脚本保留服务器上的真实项目路径和样品命名，目的是留下当前分析痕迹；换项目时直接复制对应目录并修改开头的路径、样品和contrast即可。

## 当前边界

- MACS3 peak来自每个library独立calling，再按RNA链合并；没有要求两个重复都必须call到。
- union region至少需要一个condition满足pooled IP fragments ≥5且pooled Input >0才进入检验。
- R脚本把DESeq2 dispersion固定为0，因此是Poisson Wald GLM，不是常规replicate-aware negative-binomial模型。
- 这是用于继续检验和改进m⁶A鉴定策略的`v0.1`，不是最终标准化workflow。
- 仓库不包含FASTQ、BAM、peak结果、RDS或notebook。

各项目的执行顺序见其目录中的README。
