#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsubread)
  library(DESeq2)
})

root <- Sys.getenv("ROOT", "/data/home/pengyf/projects/m6A/MTA_drought_m6A_GSE243365")
base <- file.path(root, "data", "yuhaolab_MAPQ30_PRIMARY")
peak_dir <- file.path(base, "macs3_stranded")
bam_dir <- file.path(base, "stranded")
gene_bed <- "/data/home/pengyf/reference/Araport11_TAIR10_current/genes.representative_protein_coding.bed12"
out_dir <- Sys.getenv("OUT_DIR", file.path(root, "results", "stranded_m6A_union_poisson_compact"))
threads <- as.integer(Sys.getenv("THREADS", "8"))

if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE)) > 0L &&
    Sys.getenv("OVERWRITE", "0") != "1") {
  stop("Output directory is not empty: ", out_dir)
}

say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), ..., "\n", sep="")
write_tsv <- function(x, name) write.table(
  x, file=file.path(out_dir, name), sep="\t", quote=FALSE,
  row.names=FALSE, col.names=TRUE, na="NA"
)

groups <- c("WT_Unstressed", "WT_Drought", "mta_Unstressed", "mta_Drought")
reps <- paste0("R", 1:3)
directions <- c("fwd", "rev")
design <- expand.grid(
  group=groups, replicate=reps, assay=c("Input", "m6A_IP"),
  KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE
)
design$sample_id <- paste(design$group, design$replicate, sep="_")
design$library <- paste(design$group, design$assay, design$replicate, sep="_")
stopifnot(nrow(design) == 24L, !anyDuplicated(design$library))

bam_paths <- function(direction) file.path(
  bam_dir, paste0(design$library, ".star.unique.20k.", direction, ".bam")
)
for (direction in directions) {
  missing <- bam_paths(direction)[!file.exists(bam_paths(direction))]
  if (length(missing)) stop("Missing BAMs:\n", paste(missing, collapse="\n"))
}

read_peaks <- function(path, direction) {
  if (!file.exists(path)) stop("Missing peak file: ", path)
  if (file.size(path) == 0L) return(GRanges())
  x <- read.delim(path, header=FALSE, sep="\t", quote="", comment.char="")
  GRanges(
    seqnames=x[[1]], ranges=IRanges(x[[2]] + 1L, x[[3]]),
    strand=if (direction == "fwd") "+" else "-"
  )
}

peaks <- list()
for (group in groups) for (replicate in reps) for (direction in directions) {
  key <- paste0(group, "_", replicate, ".", direction)
  path <- file.path(peak_dir, key, paste0(key, "_peaks.narrowPeak"))
  peaks[[key]] <- read_peaks(path, direction)
}
say("Loaded 24 MACS3 peak sets; raw peaks = ", sum(lengths(peaks)))

reference_union <- sort(reduce(Reduce(c, peaks), ignore.strand=FALSE))
names(reference_union) <- sprintf("union_%06d", seq_along(reference_union))
say("Union peaks = ", length(reference_union))

bed <- read.delim(gene_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
gene_ids_bed <- sub("\\.[0-9]+$", "", bed[[4]])
exon_list <- lapply(seq_len(nrow(bed)), function(i) {
  sizes <- as.integer(strsplit(as.character(bed[[11]][i]), ",", fixed=TRUE)[[1]])
  offsets <- as.integer(strsplit(as.character(bed[[12]][i]), ",", fixed=TRUE)[[1]])
  starts <- as.integer(bed[[2]][i]) + offsets + 1L
  data.frame(
    chr=bed[[1]][i], start=starts, end=starts + sizes - 1L,
    strand=bed[[6]][i], gene_id=gene_ids_bed[i]
  )
})
exon_df <- do.call(rbind, exon_list)
gene_exons <- GRanges(
  exon_df$chr, IRanges(exon_df$start, exon_df$end),
  strand=exon_df$strand, gene_id=exon_df$gene_id
)
gene_ids <- sort(unique(exon_df$gene_id))
gene_counts <- matrix(
  0L, nrow=length(gene_ids), ncol=nrow(design),
  dimnames=list(gene_ids, design$library)
)

for (direction in directions) {
  rna_strand <- if (direction == "fwd") "+" else "-"
  exons <- gene_exons[as.character(strand(gene_exons)) == rna_strand]
  saf <- data.frame(
    GeneID=mcols(exons)$gene_id, Chr=seqnames(exons),
    Start=start(exons), End=end(exons), Strand=strand(exons)
  )
  say("featureCounts gene exons: ", direction)
  fc <- featureCounts(
    files=bam_paths(direction), annot.ext=saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=TRUE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, strandSpecific=0,
    countMultiMappingReads=FALSE, ignoreDup=FALSE,
    allowMultiOverlap=FALSE, nthreads=threads
  )
  gene_counts[rownames(fc$counts), ] <- fc$counts
}

gene_size_factors <- setNames(numeric(nrow(design)), design$library)
for (assay in c("Input", "m6A_IP")) {
  libraries <- design$library[design$assay == assay]
  sf <- estimateSizeFactorsForMatrix(gene_counts[, libraries, drop=FALSE], type="poscounts")
  gene_size_factors[libraries] <- sf / exp(mean(log(sf)))
}
say("Estimated gene-exon size factors separately for Input and IP")

# Mean normalized IP/Input and log2 enrichment for each condition.
# The same calculation is applied to whole-gene exon counts and union-peak counts.
m6a_level_table <- function(counts) {
  normalized <- sweep(counts, 2, gene_size_factors[colnames(counts)], "/")
  answer <- data.frame(row.names=rownames(counts))
  for (group in groups) {
    input <- rowMeans(normalized[, paste0(group, "_Input_", reps), drop=FALSE])
    ip <- rowMeans(normalized[, paste0(group, "_m6A_IP_", reps), drop=FALSE])
    answer[[paste0(group, "_normalized_Input_mean")]] <- input
    answer[[paste0(group, "_normalized_m6A_IP_mean")]] <- ip
    answer[[paste0(group, "_log2_IP_over_Input")]] <- log2((ip + 0.5) / (input + 0.5))
  }
  answer
}

gene_exon_length <- tapply(width(gene_exons), mcols(gene_exons)$gene_id, sum)
gene_levels <- data.frame(
  gene_id=rownames(gene_counts),
  representative_exon_length=as.integer(gene_exon_length[rownames(gene_counts)]),
  m6a_level_table(gene_counts), row.names=NULL, check.names=FALSE
)

union_counts <- matrix(
  0L, nrow=length(reference_union), ncol=nrow(design),
  dimnames=list(names(reference_union), design$library)
)
union_fc_stats <- list()
for (direction in directions) {
  rna_strand <- if (direction == "fwd") "+" else "-"
  idx <- which(as.character(strand(reference_union)) == rna_strand)
  gr <- reference_union[idx]
  saf <- data.frame(
    GeneID=names(gr), Chr=seqnames(gr), Start=start(gr),
    End=end(gr), Strand=strand(gr)
  )
  say("featureCounts union peaks: ", direction, "; n = ", length(gr))
  fc <- featureCounts(
    files=bam_paths(direction), annot.ext=saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=FALSE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, strandSpecific=0,
    countMultiMappingReads=FALSE, ignoreDup=FALSE,
    allowMultiOverlap=FALSE, nthreads=threads
  )
  union_counts[idx, ] <- fc$counts[match(names(gr), rownames(fc$counts)), , drop=FALSE]
  union_fc_stats[[direction]] <- fc$stat
}

union_support <- matrix(
  FALSE, nrow=nrow(union_counts), ncol=length(groups),
  dimnames=list(rownames(union_counts), groups)
)
support_qc <- do.call(rbind, lapply(groups, function(group) {
  ip <- rowSums(union_counts[, paste0(group, "_m6A_IP_", reps), drop=FALSE])
  input <- rowSums(union_counts[, paste0(group, "_Input_", reps), drop=FALSE])
  pass <- ip >= 5 & input > 0
  union_support[, group] <<- pass
  data.frame(group=group, total=nrow(union_counts), pooled_pass=sum(pass))
}))
union_pass <- rowSums(union_support) > 0
counts_tested <- union_counts[union_pass, , drop=FALSE]
say("Union supported = ", sum(union_pass), "; unsupported = ", sum(!union_pass))

meta <- design[match(colnames(counts_tested), design$library), ]
rownames(meta) <- meta$library
meta$sample_id <- factor(meta$sample_id)
meta$assay <- factor(meta$assay, levels=c("Input", "m6A_IP"))
meta$WT_Drought_IP <- as.integer(meta$group == "WT_Drought" & meta$assay == "m6A_IP")
meta$mta_Unstressed_IP <- as.integer(meta$group == "mta_Unstressed" & meta$assay == "m6A_IP")
meta$mta_Drought_IP <- as.integer(meta$group == "mta_Drought" & meta$assay == "m6A_IP")
formula_u <- ~ sample_id + assay + WT_Drought_IP + mta_Unstressed_IP + mta_Drought_IP
X <- model.matrix(formula_u, meta)
if (qr(X)$rank != ncol(X)) stop("Design matrix is not full rank")

dds <- DESeqDataSetFromMatrix(counts_tested, meta, formula_u)
sizeFactors(dds) <- gene_size_factors[colnames(counts_tested)]
dispersions(dds) <- rep(0, nrow(dds))
say("Fitting Poisson Wald GLM: ", nrow(dds), " peaks x ", ncol(dds), " libraries")
dds <- nbinomWaldTest(dds, betaPrior=FALSE, quiet=TRUE)

native_result <- function(numerator, denominator=NULL) {
  if (is.null(denominator)) {
    result <- results(dds, name=numerator, independentFiltering=FALSE, cooksCutoff=FALSE)
  } else {
    result <- results(dds, contrast=list(numerator, denominator),
                      independentFiltering=FALSE, cooksCutoff=FALSE)
  }
  as.data.frame(result)
}

native <- list(
  WT_Drought_vs_WT_Unstressed=native_result("WT_Drought_IP"),
  mta_Drought_vs_mta_Unstressed=native_result("mta_Drought_IP", "mta_Unstressed_IP"),
  mta_Unstressed_vs_WT_Unstressed=native_result("mta_Unstressed_IP"),
  mta_Drought_vs_WT_Drought=native_result("mta_Drought_IP", "WT_Drought_IP")
)

master <- data.frame(
  peak_id=names(reference_union), chr=seqnames(reference_union),
  start=start(reference_union), end=end(reference_union),
  strand=strand(reference_union), width_bp=width(reference_union),
  read_supported=union_pass, stringsAsFactors=FALSE
)
for (group in groups) master[[paste0("support_", group)]] <- union_support[, group]
for (library in design$library) master[[paste0("count_", library)]] <- union_counts[, library]
peak_levels <- m6a_level_table(union_counts)
master <- cbind(master, peak_levels[master$peak_id, , drop=FALSE])

paper_call <- function(lfc, pvalue, tested, fc=1.5) {
  answer <- ifelse(tested, "Not_selected", "Not_tested")
  answer[tested & is.finite(lfc) & is.finite(pvalue) & pvalue < .05 & lfc >= log2(fc)] <- "Up"
  answer[tested & is.finite(lfc) & is.finite(pvalue) & pvalue < .05 & lfc <= -log2(fc)] <- "Down"
  answer
}

for (contrast in names(native)) {
  result <- native[[contrast]]
  idx <- match(master$peak_id, rownames(result))
  for (column in colnames(result)) master[[paste0(contrast, "_", column)]] <- result[[column]][idx]
  master[[paste0(contrast, "_call_FC1.5_P0.05")]] <- paper_call(
    master[[paste0(contrast, "_log2FoldChange")]],
    master[[paste0(contrast, "_pvalue")]], master$read_supported
  )
}

summary_rows <- list()
for (fc in c(1.5, 2)) for (contrast in names(native)) {
  lfc <- master[[paste0(contrast, "_log2FoldChange")]]
  pvalue <- master[[paste0(contrast, "_pvalue")]]
  selected <- is.finite(lfc) & is.finite(pvalue) & pvalue < .05
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    contrast=contrast, FC_threshold=fc, P_threshold=.05,
    Up=sum(selected & lfc >= log2(fc)),
    Down=sum(selected & lfc <= -log2(fc))
  )
}
summary_tab <- do.call(rbind, summary_rows)
print(summary_tab, row.names=FALSE)

library_qc <- data.frame(
  library=design$library, group=design$group, replicate=design$replicate,
  assay=design$assay, gene_size_factor=gene_size_factors[design$library],
  assigned_gene_fragments=colSums(gene_counts),
  assigned_union_fragments=colSums(union_counts),
  union_peaks_with_reads=colSums(union_counts > 0)
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write_tsv(master, "union_poisson_all_peaks.tsv")
write_tsv(gene_levels, "gene_exon_m6A_levels.tsv")
write.table(
  data.frame(chr=master$chr, start0=master$start - 1L, end=master$end,
             peak_id=master$peak_id, score=0L, strand=master$strand),
  file=file.path(out_dir, "union_peaks.bed"), sep="\t", quote=FALSE,
  row.names=FALSE, col.names=FALSE
)
write_tsv(library_qc, "library_qc_size_factors.tsv")
write_tsv(support_qc, "read_support_qc.tsv")
write_tsv(summary_tab, "differential_summary.tsv")

saveRDS(list(
  reference_union=reference_union, union_counts=union_counts,
  union_support=union_support, union_pass=union_pass,
  gene_counts=gene_counts, gene_size_factors=gene_size_factors,
  gene_levels=gene_levels, peak_levels=peak_levels,
  design=design, meta=meta, model_matrix=X, dds=dds,
  native_results=native, master=master, summary=summary_tab,
  session_info=sessionInfo()
), file.path(out_dir, "union_poisson_objects.rds"), compress="xz")

writeLines(c(
  "GSE243365 | stranded MACS3 union + paired Input/IP Poisson GLM",
  paste("Date:", format(Sys.time(), "%F %T %Z")),
  "Inputs: MAPQ30 primary, 20-kb TLEN filter, reverse-strand split BAMs.",
  "Union: 24 sample/strand MACS3 peak sets; fwd/rev remain separate.",
  "Counts: paired fragments; gene-exon size factors estimated separately for Input and IP.",
  "m6A level: log2((mean normalized IP + 0.5)/(mean normalized Input + 0.5)).",
  "Peak-level values are columns in union_poisson_all_peaks.tsv.",
  "Whole-gene representative-exon values are in gene_exon_m6A_levels.tsv.",
  "Support: pooled three-replicate IP >=5 and pooled Input >0 in at least one group.",
  "Model: ~ sample_id + assay + WT_Drought_IP + mta_Unstressed_IP + mta_Drought_IP.",
  "Dispersion is fixed to zero to reproduce the exploratory Poisson branch.",
  "Selection columns use raw P<0.05 and |FC|>=1.5; native BH padj is also retained.",
  "This is not a complete exomePeak2 analysis and Poisson significance can be anti-conservative."
), file.path(out_dir, "README.txt"))

say("Saved: ", out_dir)
