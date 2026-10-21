#!/usr/bin/env Rscript
# PRJNA1028096: stranded MACS3 union peaks followed by paired Poisson GLMs.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsubread)
  library(DESeq2)
})

root <- "/data/home/pengyf/projects/m6A/Defining_context_dependent_m6A_PRJNA1028096"
base <- file.path(root, "data", "yuhaolab_MAPQ30_PRIMARY_SECOND_STRAND")
sample_file <- file.path(root, "config", "sample_exomepeak2_temperature_and_light_quality_mapq30_primary.txt")
contrast_file <- file.path(root, "config", "compare_exomepeak2_temperature_and_light_quality.tsv")
annotation_file <- "/data/home/pengyf/reference/Araport11_TAIR10_current/genes.representative_protein_coding.bed12"
out_dir <- file.path(root, "results", "stranded_m6A_union_poisson_compact")
threads <- as.integer(Sys.getenv("THREADS", "8"))

samples <- read.delim(sample_file, check.names=FALSE)
contrasts <- read.delim(contrast_file, check.names=FALSE)
samples$assay <- ifelse(grepl("_IP_", samples$sample), "IP", "Input")
samples$replicate <- sub(".*_(R[0-9]+)$", "\\1", samples$sample)
samples$library <- samples$sample
stopifnot(nrow(samples) == 52L, !anyDuplicated(samples$library))

bam_dir <- file.path(base, "stranded")
peak_dir <- file.path(base, "macs3_stranded")
directions <- c("fwd", "rev")
bam_paths <- function(direction, libraries=samples$library) {
  file.path(bam_dir, paste0(libraries, ".star.unique.20k.", direction, ".bam"))
}

read_peaks <- function(sample, direction) {
  name <- paste(sample, direction, sep=".")
  path <- file.path(peak_dir, name, paste0(name, "_peaks.narrowPeak"))
  if (!file.exists(path)) stop("Missing MACS3 peaks: ", path)
  if (file.size(path) == 0L) return(GRanges())
  x <- read.delim(path, header=FALSE, sep="\t")
  GRanges(x[[1]], IRanges(x[[2]] + 1L, x[[3]]),
          strand=if (direction == "fwd") "+" else "-")
}

ip_samples <- samples$sample[samples$assay == "IP"]
peaks <- unlist(GRangesList(lapply(ip_samples, function(s) {
  c(read_peaks(s, "fwd"), read_peaks(s, "rev"))
})), use.names=FALSE)
union <- sort(reduce(peaks, ignore.strand=FALSE))
names(union) <- sprintf("union_%06d", seq_along(union))
cat("Union peaks:", length(union), "\n")

# Representative protein-coding exons for robust library-size normalization.
bed <- read.delim(annotation_file, header=FALSE, sep="\t")
exons <- do.call(rbind, lapply(seq_len(nrow(bed)), function(i) {
  size <- as.integer(strsplit(as.character(bed[i, 11]), ",", fixed=TRUE)[[1]])
  offset <- as.integer(strsplit(as.character(bed[i, 12]), ",", fixed=TRUE)[[1]])
  start <- bed[i, 2] + offset + 1L
  data.frame(gene_id=sub("\\.[0-9]+$", "", bed[i, 4]), chr=bed[i, 1],
             start=start, end=start + size - 1L, strand=bed[i, 6])
}))
gene_ids <- sort(unique(exons$gene_id))
gene_counts <- matrix(0L, length(gene_ids), nrow(samples),
                      dimnames=list(gene_ids, samples$library))
union_counts <- matrix(0L, length(union), nrow(samples),
                       dimnames=list(names(union), samples$library))

for (direction in directions) {
  strand_value <- if (direction == "fwd") "+" else "-"
  exon_saf <- transform(exons[exons$strand == strand_value, ],
                        GeneID=gene_id, Chr=chr, Start=start, End=end, Strand=strand_value)
  exon_saf <- exon_saf[, c("GeneID", "Chr", "Start", "End", "Strand")]
  fc_gene <- featureCounts(
    files=bam_paths(direction), annot.ext=exon_saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=TRUE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, countMultiMappingReads=FALSE,
    allowMultiOverlap=FALSE, nthreads=threads
  )
  gene_counts[rownames(fc_gene$counts), ] <- fc_gene$counts

  idx <- which(as.character(strand(union)) == strand_value)
  gr <- union[idx]
  peak_saf <- data.frame(GeneID=names(gr), Chr=as.character(seqnames(gr)),
                         Start=start(gr), End=end(gr), Strand=strand_value)
  fc_peak <- featureCounts(
    files=bam_paths(direction), annot.ext=peak_saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=FALSE,
    isPairedEnd=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE,
    countChimericFragments=FALSE, countMultiMappingReads=FALSE,
    allowMultiOverlap=FALSE, nthreads=threads
  )
  union_counts[idx, ] <- fc_peak$counts[match(names(gr), rownames(fc_peak$counts)), ]
}

# A union region is tested if at least one condition has pooled IP >=5 and Input >0.
groups <- unique(samples$group)
support <- sapply(groups, function(group) {
  ip <- samples$library[samples$group == group & samples$assay == "IP"]
  input <- samples$library[samples$group == group & samples$assay == "Input"]
  rowSums(union_counts[, ip, drop=FALSE]) >= 5 &
    rowSums(union_counts[, input, drop=FALSE]) > 0
})
colnames(support) <- groups
tested <- rowSums(support) > 0

master <- data.frame(
  peak_id=names(union), chr=as.character(seqnames(union)),
  start=start(union), end=end(union), strand=as.character(strand(union)),
  read_supported=tested, stringsAsFactors=FALSE
)
for (library in samples$library) master[[paste0("count_", library)]] <- union_counts[, library]

fit_contrast <- function(reference, treatment) {
  keep_samples <- samples$group %in% c(reference, treatment)
  meta <- samples[keep_samples, ]
  counts <- union_counts[tested, meta$library, drop=FALSE]
  keep_peaks <- rowSums(counts) > 0
  counts <- counts[keep_peaks, , drop=FALSE]
  meta$pair <- factor(paste(meta$group, meta$replicate, sep="_"))
  meta$assay <- factor(meta$assay, levels=c("Input", "IP"))
  meta$treatment_IP <- as.integer(meta$group == treatment & meta$assay == "IP")
  rownames(meta) <- meta$library

  size_factor <- numeric(nrow(meta))
  names(size_factor) <- meta$library
  for (assay in c("Input", "IP")) {
    libraries <- meta$library[meta$assay == assay]
    value <- estimateSizeFactorsForMatrix(
      gene_counts[, libraries, drop=FALSE], type="poscounts"
    )
    size_factor[libraries] <- value / exp(mean(log(value)))
  }

  dds <- DESeqDataSetFromMatrix(counts, meta, ~ pair + assay + treatment_IP)
  sizeFactors(dds) <- size_factor[colnames(dds)]
  dispersions(dds) <- 0
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, quiet=TRUE)
  result <- as.data.frame(results(dds, name="treatment_IP",
                                  independentFiltering=FALSE,
                                  cooksCutoff=FALSE))
  result$peak_id <- rownames(result)
  result
}

all_results <- list()
for (i in seq_len(nrow(contrasts))) {
  reference <- contrasts$group1[i]
  treatment <- contrasts$group2[i]
  name <- paste(treatment, "vs", reference, sep="_")
  cat("Fitting", name, "\n")
  result <- fit_contrast(reference, treatment)
  all_results[[name]] <- result
  idx <- match(master$peak_id, result$peak_id)
  for (column in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    master[[paste0(name, "_", column)]] <- result[[column]][idx]
  }
}

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write.table(master, file.path(out_dir, "union_poisson_all_peaks.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(data.frame(chr=master$chr, start=master$start - 1L,
                       end=master$end, peak_id=master$peak_id,
                       score=0, strand=master$strand),
            file.path(out_dir, "union_peaks.bed"), sep="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
saveRDS(list(union=union, counts=union_counts, support=support,
             gene_counts=gene_counts, samples=samples, contrasts=contrasts,
             results=all_results, master=master, session_info=sessionInfo()),
        file.path(out_dir, "union_poisson_objects.rds"), compress="xz")
cat("Wrote results to", out_dir, "\n")
