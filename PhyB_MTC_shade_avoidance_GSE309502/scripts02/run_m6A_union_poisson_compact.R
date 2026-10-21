#!/usr/bin/env Rscript
# GSE309502 | Li et al. 2026 Figure 4-style m6A reanalysis
# Final path taken from run02.stranded_m6A_QC_and_counts(3).ipynb:
#   MACS3 fwd/rev union -> gene-level size factors -> stranded peak counts
#   -> support filter -> paired DESeq2 Poisson Wald GLM
# Not a direct exomePeak2 invocation; not a complete exomePeak2 workflow.
#
# Outputs: ONE all-union peak table, ONE comprehensive RDS object, plus
#          union BED, library QC, support QC, brief FC/P summary, README.
# No exploratory NB, comparisons/plots, gene annotation, or candidate rankings.
# Original BAMs and peak sets are never modified.
#
# Run from the project root:
#   THREADS=8 Rscript run_m6A_union_poisson_compact.R > run_union_poisson.log 2>&1
# Use OUT_DIR=... to change destination; by default nonempty outputs cause stop.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsubread)
  library(DESeq2)
})

root <- Sys.getenv(
  "ROOT", "/data/home/pengyf/projects/m6A/PhyB_MTC_shade_avoidance_GSE309502"
)
base <- file.path(root, "data", "yuhaolab_MAPQ30_PRIMARY")
peak_dir <- file.path(base, "macs3_stranded")
bam_dir <- file.path(base, "stranded")
annotation_file <- Sys.getenv(
  "GENE_BED12",
  "/data/home/pengyf/reference/Araport11_TAIR10_current/genes.representative_protein_coding.bed12"
)
out_dir <- Sys.getenv("OUT_DIR", file.path(root, "results", "stranded_m6A_union_poisson_compact"))
threads <- as.integer(Sys.getenv("THREADS", "8"))
if (is.na(threads) || threads < 1L) stop("THREADS must be >= 1")
if (!file.exists(annotation_file)) stop("Missing Araport11 BED12: ", annotation_file)
if (!dir.exists(peak_dir)) stop("Missing MACS3 peak directory: ", peak_dir)
if (!dir.exists(bam_dir)) stop("Missing stranded BAM directory: ", bam_dir)

# Do not overwrite an existing set of results by accident.
if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE)) > 0L &&
    Sys.getenv("OVERWRITE", "0") != "1") {
  stop("Output directory is not empty: ", out_dir,
       "\nSpecify a new OUT_DIR, or set OVERWRITE=1 to replace matching files.")
}

say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), ..., "\n", sep = "")
write_tsv <- function(x, filename) {
  write.table(x, file=file.path(out_dir, filename), sep="\t", quote=FALSE,
              row.names=FALSE, col.names=TRUE, na="NA")
}

# -----------------------------------------------------------------------------
# 1. Three conditions x two biological replicates x Input/IP = 12 libraries
# -----------------------------------------------------------------------------
groups <- c("Col-0_WL", "Col-0_Shade", "phyB-9_WL")
reps <- c("R1", "R2")
directions <- c("fwd", "rev")
design <- expand.grid(
  group=groups, replicate=reps, assay=c("Input", "IP"),
  KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE
)
design$sample_id <- paste(design$group, design$replicate, sep="_")
design$library <- paste(design$group, design$assay, design$replicate, sep="_")
stopifnot(nrow(design) == 12L, !anyDuplicated(design$library))

bam_paths <- function(direction) file.path(
  bam_dir, paste0(design$library, ".star.unique.20k.", direction, ".bam")
)
for (direction in directions) {
  bams <- bam_paths(direction)
  if (!all(file.exists(bams)))
    stop("Missing ", direction, " BAM:\n", paste(bams[!file.exists(bams)], collapse="\n"))
}

# -----------------------------------------------------------------------------
# 2. Read 12 MACS3 narrowPeak sets and merge overlaps on the SAME RNA strand
#    BED: 0-based half-open --> GRanges: 1-based closed
# -----------------------------------------------------------------------------
read_peaks <- function(path, direction) {
  if (!file.exists(path)) stop("Missing peak file: ", path)
  if (file.size(path) == 0L) return(GRanges())
  x <- read.delim(path, header=FALSE, sep="\t", quote="",
                  comment.char="", stringsAsFactors=FALSE)
  if (ncol(x) < 9L) stop("Invalid narrowPeak file: ", path)
  if (any(x[[2]] < 0 | x[[3]] <= x[[2]])) stop("Invalid peak coordinates: ", path)
  GRanges(
    seqnames=x[[1]],
    ranges=IRanges(start=as.integer(x[[2]]) + 1L, end=as.integer(x[[3]])),
    strand=if (direction == "fwd") "+" else "-",
    signalValue=as.numeric(x[[7]]), negLog10P=as.numeric(x[[8]])
  )
}
peaks <- list()
for (g in groups) for (r in reps) for (d in directions) {
  key <- paste0(g, "_", r, ".", d)
  path <- file.path(peak_dir, key, paste0(key, "_peaks.narrowPeak"))
  peaks[[key]] <- read_peaks(path, d)
}
say("Loaded 12 MACS3 peak sets, raw peaks total = ", sum(lengths(peaks)))

reference_union <- sort(reduce(Reduce(c, peaks), ignore.strand=FALSE))
names(reference_union) <- sprintf("union_%06d", seq_along(reference_union))
stopifnot(!anyDuplicated(names(reference_union)),
          all(as.character(strand(reference_union)) %in% c("+", "-")))
say("Union peaks = ", length(reference_union),
    "; + ", sum(strand(reference_union) == "+"),
    "; - ", sum(strand(reference_union) == "-"))

# -----------------------------------------------------------------------------
# 3. Count fragments on gene exons to estimate library size factors.
#    Matches 02C.1 + 02E.1, with separate IP and Input normalization.
# -----------------------------------------------------------------------------
bed <- read.delim(annotation_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
if (ncol(bed) < 12L) stop("BED12 must contain at least 12 columns")
gene_ids_bed <- sub("\\.[0-9]+$", "", bed[[4]])

exon_list <- lapply(seq_len(nrow(bed)), function(i) {
  sizes <- as.integer(strsplit(as.character(bed[[11]][i]), ",", fixed=TRUE)[[1]])
  offsets <- as.integer(strsplit(as.character(bed[[12]][i]), ",", fixed=TRUE)[[1]])
  if (length(sizes) != length(offsets) || anyNA(sizes) || anyNA(offsets))
    stop("Invalid exon blocks at annotation row ", i)
  exon_start <- as.integer(bed[[2]][i]) + offsets + 1L
  data.frame(chr=bed[[1]][i], start=exon_start,
             end=exon_start + sizes - 1L,
             strand=bed[[6]][i], gene_id=gene_ids_bed[i],
             stringsAsFactors=FALSE)
})
exon_df <- do.call(rbind, exon_list)
gene_exons <- GRanges(seqnames=exon_df$chr,
                      ranges=IRanges(exon_df$start, exon_df$end),
                      strand=exon_df$strand, gene_id=exon_df$gene_id)
gene_ids <- sort(unique(as.character(mcols(gene_exons)$gene_id)))
gene_counts <- matrix(0L, nrow=length(gene_ids), ncol=nrow(design),
                      dimnames=list(gene_ids, design$library))

for (d in directions) {
  s <- if (d == "fwd") "+" else "-"
  exons <- gene_exons[as.character(strand(gene_exons)) == s]
  saf <- data.frame(
    GeneID=as.character(mcols(exons)$gene_id),
    Chr=as.character(seqnames(exons)),
    Start=start(exons), End=end(exons), Strand=as.character(strand(exons))
  )
  say("featureCounts gene exons, ", d, ", entries = ", nrow(saf))
  fc <- Rsubread::featureCounts(
    files=bam_paths(d), annot.ext=saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=TRUE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, strandSpecific=0,
    countMultiMappingReads=FALSE, ignoreDup=FALSE, allowMultiOverlap=FALSE,
    nthreads=threads
  )
  stopifnot(!anyDuplicated(rownames(fc$counts)),
            all(rownames(fc$counts) %in% gene_ids))
  gene_counts[rownames(fc$counts), ] <- fc$counts
}

gene_size_factors <- setNames(numeric(nrow(design)), design$library)
for (a in c("Input", "IP")) {
  samples <- design$library[design$assay == a]
  z <- DESeq2::estimateSizeFactorsForMatrix(
    gene_counts[, samples, drop=FALSE], type="poscounts"
  )
  z <- z / exp(mean(log(z))) # geom. mean = 1 separately within each assay
  gene_size_factors[samples] <- z
}
stopifnot(all(is.finite(gene_size_factors)), all(gene_size_factors > 0))
say("Estimated gene-level size factors separately for Input and IP")

# -----------------------------------------------------------------------------
# 4. Count fwd union peaks against fwd BAMs, rev against rev BAMs.
#    One paired-end fragment per count; identical featureCounts settings as 02F.
# -----------------------------------------------------------------------------
union_counts <- matrix(0L, nrow=length(reference_union), ncol=nrow(design),
                       dimnames=list(names(reference_union), design$library))
union_fc_stats <- list()
for (d in directions) {
  s <- if (d == "fwd") "+" else "-"
  idx <- which(as.character(strand(reference_union)) == s)
  gr <- reference_union[idx]
  saf <- data.frame(GeneID=names(gr), Chr=as.character(seqnames(gr)),
                    Start=start(gr), End=end(gr), Strand=as.character(strand(gr)))
  say("featureCounts union peaks, ", d, ", peaks = ", length(idx))
  fc <- Rsubread::featureCounts(
    files=bam_paths(d), annot.ext=saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=FALSE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, strandSpecific=0,
    countMultiMappingReads=FALSE, ignoreDup=FALSE, allowMultiOverlap=FALSE,
    nthreads=threads
  )
  ix <- match(saf$GeneID, rownames(fc$counts))
  if (anyNA(ix)) stop("Union featureCounts IDs did not match ", d)
  union_counts[idx, ] <- fc$counts[ix, , drop=FALSE]
  union_fc_stats[[d]] <- fc$stat
}
stopifnot(identical(rownames(union_counts), names(reference_union)), any(union_counts > 0))

# Identical support rule to 02F.2: pooled IP>=5 AND pooled Input>0 in >=1 group.
union_support <- matrix(FALSE, nrow=nrow(union_counts), ncol=length(groups),
                        dimnames=list(rownames(union_counts), groups))
union_filter_qc <- do.call(rbind, lapply(groups, function(g) {
  ip1 <- union_counts[, paste0(g, "_IP_R1")]
  ip2 <- union_counts[, paste0(g, "_IP_R2")]
  in1 <- union_counts[, paste0(g, "_Input_R1")]
  in2 <- union_counts[, paste0(g, "_Input_R2")]
  pooled <- (ip1 + ip2 >= 5) & (in1 + in2 > 0)
  r1 <- (ip1 >= 5) & (in1 > 0)
  r2 <- (ip2 >= 5) & (in2 > 0)
  union_support[, g] <<- pooled
  data.frame(group=g, total=nrow(union_counts), pooled_pass=sum(pooled),
             both_reps_pass=sum(r1 & r2))
}))
rownames(union_filter_qc) <- NULL
union_pass <- rowSums(union_support) > 0
union_reference_supported <- reference_union[union_pass]
union_counts_supported <- union_counts[union_pass, , drop=FALSE]
stopifnot(nrow(union_counts_supported) > 0L)
say("Union supported = ", sum(union_pass), "; not supported = ", sum(!union_pass))

# -----------------------------------------------------------------------------
# 5. DESeq2 Poisson (NB GLM with dispersion fixed to zero).
#    This reproduces the notebook's Poisson BRANCH, not exomePeak2 end-to-end.
# -----------------------------------------------------------------------------
meta <- design[match(colnames(union_counts_supported), design$library), , drop=FALSE]
stopifnot(!anyNA(meta$library), identical(meta$library, colnames(union_counts_supported)))
rownames(meta) <- meta$library
meta$sample_id <- factor(meta$sample_id)
meta$assay <- factor(meta$assay, levels=c("Input", "IP"))
meta$shade_IP <- as.integer(meta$group == "Col-0_Shade" & meta$assay == "IP")
meta$phyB_IP <- as.integer(meta$group == "phyB-9_WL" & meta$assay == "IP")
formula_u <- ~ sample_id + assay + shade_IP + phyB_IP
X <- model.matrix(formula_u, data=meta)
if (ncol(X) != 9L || qr(X)$rank != ncol(X)) stop("Design matrix is not full rank")

dds_pois_union <- DESeq2::DESeqDataSetFromMatrix(
  countData=union_counts_supported, colData=meta, design=formula_u
)
sf <- gene_size_factors[colnames(union_counts_supported)]
stopifnot(all(is.finite(sf)), all(sf > 0))
sizeFactors(dds_pois_union) <- as.numeric(sf)
dispersions(dds_pois_union) <- rep(0, nrow(dds_pois_union))
say("Fitting Poisson Wald test, 21k peaks x 12 libraries")
dds_pois_union <- DESeq2::nbinomWaldTest(dds_pois_union, betaPrior=FALSE, quiet=TRUE)
stopifnot(all(c("shade_IP", "phyB_IP") %in% resultsNames(dds_pois_union)))

get_native <- function(coef_name) {
  ans <- DESeq2::results(dds_pois_union, name=coef_name,
                         independentFiltering=FALSE, cooksCutoff=FALSE)
  ans <- as.data.frame(ans)  # baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
  stopifnot(identical(rownames(ans), rownames(union_counts_supported)))
  ans
}
res_shade <- get_native("shade_IP")
res_phyB <- get_native("phyB_IP")

# -----------------------------------------------------------------------------
# 6. Restore notebook-compatible result structures; build ONE complete table.
#    Raw DESeq2 output is retained in BOTH the table and RDS.
# -----------------------------------------------------------------------------
union_results <- list(
  Poisson = list(
    Shade = data.frame(log2FC=res_shade$log2FoldChange,
                       pvalue=res_shade$pvalue,
                       fdr=res_shade$padj,
                       row.names=rownames(res_shade)),
    phyB = data.frame(log2FC=res_phyB$log2FoldChange,
                      pvalue=res_phyB$pvalue,
                      fdr=res_phyB$padj,
                      row.names=rownames(res_phyB))
  )
)
poisson_native_results <- list(Shade=res_shade, phyB=res_phyB)

# 21,797 rows including unsupported union regions. Unsupported statistics = NA.
master <- data.frame(
  peak_id=names(reference_union),
  chr=as.character(seqnames(reference_union)),
  start=start(reference_union), end=end(reference_union),
  strand=as.character(strand(reference_union)),
  width_bp=width(reference_union),
  read_supported=as.logical(union_pass),
  stringsAsFactors=FALSE
)
for (g in groups) master[[paste0("support_", g)]] <- as.logical(union_support[, g])
for (lib in design$library) master[[paste0("count_", lib)]] <- union_counts[, lib]

for (contrast in c("Shade", "phyB")) {
  res <- poisson_native_results[[contrast]]
  idx <- match(master$peak_id, rownames(res))
  for (nm in colnames(res)) {
    master[[paste0(contrast, "_", nm)]] <- res[[nm]][idx]
  }
}

paper_call <- function(lfc, p, tested, fc=1.5, p_cutoff=0.05) {
  ans <- rep("Not_selected", length(lfc))
  ans[!tested] <- "Not_tested"
  ok <- tested & is.finite(lfc) & is.finite(p) & (p < p_cutoff)
  ans[ok & lfc >= log2(fc)] <- "Up"
  ans[ok & lfc <= -log2(fc)] <- "Down"
  ans
}
master$Shade_call_FC1.5_P0.05 <- paper_call(
  master$Shade_log2FoldChange, master$Shade_pvalue, master$read_supported
)
master$phyB_call_FC1.5_P0.05 <- paper_call(
  master$phyB_log2FoldChange, master$phyB_pvalue, master$read_supported
)
master$common_call_FC1.5_P0.05 <- ifelse(
  master$Shade_call_FC1.5_P0.05 == "Up" &
    master$phyB_call_FC1.5_P0.05 == "Up", "Up",
  ifelse(master$Shade_call_FC1.5_P0.05 == "Down" &
         master$phyB_call_FC1.5_P0.05 == "Down", "Down", "Not_selected")
)
master$common_call_FC1.5_P0.05[!master$read_supported] <- "Not_tested"

count_direction <- function(df, contrast, fc, direction) {
  lfc <- df[[paste0(contrast, "_log2FoldChange")]]
  p <- df[[paste0(contrast, "_pvalue")]]
  valid <- is.finite(lfc) & is.finite(p) & p < 0.05
  sum(valid & if (direction == "Up") lfc >= log2(fc) else lfc <= -log2(fc))
}
summary_tab <- do.call(rbind, lapply(c(1.5, 2), function(fc) {
  su <- count_direction(master, "Shade", fc, "Up")
  sd <- count_direction(master, "Shade", fc, "Down")
  pu <- count_direction(master, "phyB", fc, "Up")
  pd <- count_direction(master, "phyB", fc, "Down")
  shade_lfc <- master$Shade_log2FoldChange
  phyB_lfc <- master$phyB_log2FoldChange
  valid <- is.finite(shade_lfc) & is.finite(phyB_lfc) &
    is.finite(master$Shade_pvalue) & is.finite(master$phyB_pvalue) &
    master$Shade_pvalue < 0.05 & master$phyB_pvalue < 0.05
  data.frame(FC_threshold=fc, P_threshold=0.05,
             Shade_Up=su, Shade_Down=sd, phyB_Up=pu, phyB_Down=pd,
             Common_Up=sum(valid & shade_lfc >= log2(fc) & phyB_lfc >= log2(fc)),
             Common_Down=sum(valid & shade_lfc <= -log2(fc) & phyB_lfc <= -log2(fc)))
}))
rownames(summary_tab) <- NULL
say("Differential peak counts:\n", paste(capture.output(print(summary_tab, row.names=FALSE)), collapse="\n"))

library_qc <- data.frame(
  library=design$library, group=design$group,
  replicate=design$replicate, assay=design$assay,
  gene_size_factor=unname(gene_size_factors[design$library]),
  assigned_gene_fragments=as.numeric(colSums(gene_counts)),
  assigned_union_fragments=as.numeric(colSums(union_counts)),
  union_peaks_with_reads=as.numeric(colSums(union_counts > 0))
)

# -----------------------------------------------------------------------------
# 7. Export only the requested files (no separate Up/Down/contrast TSVs).
# -----------------------------------------------------------------------------
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(out_dir)) stop("Cannot create OUT_DIR: ", out_dir)

# A single full native results table, not just significant peaks.
# TSV coordinates are 1-based closed; all 12 fragment counts are unmodified.
write_tsv(master, "union_poisson_all_peaks.tsv")

# BED6 for future IGV/BEDtools use; BED start is 0-based half-open.
all_bed <- data.frame(
  chr=master$chr, start0=master$start - 1L, end=master$end,
  peak_id=master$peak_id, score=0L, strand=master$strand
)
write.table(all_bed, file=file.path(out_dir, "union_peaks.bed"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write_tsv(library_qc, "library_qc_size_factors.tsv")
write_tsv(union_filter_qc, "read_support_qc.tsv")
write_tsv(summary_tab, "differential_summary.tsv")

# One complete R object, keeping notebook variable names for downstream reuse.
# Unmodified DESeq2 native result data.frames are also available in RDS.
result_objects <- list(
  reference_union=reference_union,
  union_counts=union_counts,
  union_support=union_support,
  union_pass=union_pass,
  union_reference_supported=union_reference_supported,
  union_counts_supported=union_counts_supported,
  union_fc_stats=union_fc_stats,
  gene_exons=gene_exons,
  gene_counts=gene_counts,
  gene_size_factors=gene_size_factors,
  design=design, meta=meta, formula_u=formula_u,
  model_matrix=X,
  dds_pois_union=dds_pois_union,
  union_results=union_results,
  poisson_native_results=poisson_native_results,
  union_filter_qc=union_filter_qc,
  library_qc=library_qc,
  summary=summary_tab,
  master=master,
  session_info=sessionInfo()
)
saveRDS(result_objects, file.path(out_dir, "union_poisson_objects.rds"), compress="xz")

readme <- c(
  "GSE309502 | MACS3 stranded union + DESeq2 Poisson GLM",
  "This is the FINAL Poisson branch from run02 Notebook, NOT a full exomePeak2 run.",
  paste("Date:", format(Sys.time(), "%F %T %Z")),
  paste("Input peak sets:", peak_dir),
  paste("Input stranded BAMs:", bam_dir),
  paste("Gene annotation:", annotation_file),
  "Method: same-strand MACS3 union -> fwd/rev featureCounts paired fragments ->",
  "  gene-exon DESeq2 size factors separately for Input and IP (poscounts) ->",
  "  pool R1+R2 IP >=5 and Input >0 in ANY group ->",
  "  design = ~sample_id + assay + shade_IP + phyB_IP; dispersion=0;",
  "  nbinomWaldTest(betaPrior=FALSE); results independentFiltering=FALSE, cooksCutoff=FALSE.",
  "Main threshold: raw pvalue <0.05, log2FoldChange >= log2(1.5) for Up",
  "  or log2FoldChange <= -log2(1.5) for Down; BH padj is reported but NOT used to select.",
  "Main TSV: union_poisson_all_peaks.tsv; ALL union peaks and both contrasts.",
  "  12 raw fragment-count columns; native DESeq2 baseMean, log2FoldChange,",
  "  lfcSE, stat, pvalue, padj for BOTH contrasts; unsupported regions have NA statistics.",
  "  TSV genomic positions are 1-based closed; BED is 0-based half-open.",
  "RDS: union_poisson_objects.rds; includes full dds_pois_union and original", 
  "  notebook-compatible union_results$Poisson$Shade/$phyB and all read counts.",
  "Supporting: union_peaks.bed, library_qc_size_factors.tsv,",
  "  read_support_qc.tsv, differential_summary.tsv.",
  "No peak-to-gene assignment, no detailed candidate TSVs, no IGV ranking.",
  "Caution: Poisson dispersion=0 may underestimate variability; matches total peak",
  "  counts approximately but does not validate FDR significance or individual peaks.",
  paste("Union:", nrow(master), "Tested:", sum(master$read_supported))
)
writeLines(readme, file.path(out_dir, "README.txt"))

expected <- c(21797L, 21474L, 1582L, 2118L, 1598L, 1957L, 565L, 888L)
observed <- c(nrow(master), sum(master$read_supported),
              unlist(summary_tab[1, c("Shade_Up", "Shade_Down", "phyB_Up",
                                      "phyB_Down", "Common_Up", "Common_Down")],
                     use.names=FALSE))
if (!identical(as.integer(observed), expected)) {
  warning("Observed counts differ from notebook: ", paste(observed, collapse=", "),
          "; expected: ", paste(expected, collapse=", "))
} else {
  say("All eight reference checks MATCH notebook")
}
say("Saved compact result set in: ", out_dir)
