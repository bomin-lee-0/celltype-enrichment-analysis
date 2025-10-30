#!/usr/bin/env Rscript
# ============================================================
# Peak Annotation with ChIPseeker
# Annotates peaks to genomic features using rat genome (rn7)
# ============================================================

cat("========================================\n")
cat("Peak Annotation with ChIPseeker\n")
cat("========================================\n\n")

# Set working directory
setwd("/scratch/prj/bcn_marzi_lab/ratlas/Bomin")

# Load libraries
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
    library(org.Rn.eg.db)
    library(dplyr)
    library(ggplot2)
})

cat("Libraries loaded successfully\n\n")

# Load TxDb
txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene

# Create output directory
dir.create("annotation", showWarnings = FALSE)

# ============================================================
# Input files - MACS2 narrowPeak files
# ============================================================
files <- list(
  NeuN = "macs2_output/NeuN/NeuN_all_peaks.narrowPeak",
  Nurr = "macs2_output/Nurr/Nurr_all_peaks.narrowPeak",
  Olig = "macs2_output/Olig/Olig_all_peaks.narrowPeak",
  Neg  = "macs2_output/Neg/Neg_all_peaks.narrowPeak"
)

cat("Input files:\n")
for (name in names(files)) {
    if (file.exists(files[[name]])) {
        peak_count <- length(readLines(files[[name]]))
        cat(sprintf("  ✓ %-10s: %s (%d peaks)\n", name, basename(files[[name]]), peak_count))
    } else {
        cat(sprintf("  ✗ %-10s: %s (NOT FOUND)\n", name, basename(files[[name]])))
    }
}
cat("\n")

# ============================================================
# Annotate peaks
# ============================================================
cat("[", as.character(Sys.time()), "] Starting peak annotation...\n\n")

peakAnnoList <- lapply(files, function(peak_file) {
    if (!file.exists(peak_file)) {
        return(NULL)
    }

    cat("Annotating:", basename(peak_file), "\n")

    anno <- annotatePeak(
        peak_file,
        TxDb = txdb,
        tssRegion = c(-3000, 3000),
        annoDb = "org.Rn.eg.db"
    )

    cat("  Peaks annotated:", length(anno@anno), "\n")
    cat("  Genes identified:", length(unique(anno@anno$geneId)), "\n\n")

    return(anno)
})

# Remove NULL entries
peakAnnoList <- peakAnnoList[!sapply(peakAnnoList, is.null)]

cat("[", as.character(Sys.time()), "] Annotation complete\n\n")

# ============================================================
# Convert to data frames
# ============================================================
cat("Converting annotations to data frames...\n")

annotation_summary <- lapply(peakAnnoList, function(x) {
    if (!is.null(x)) {
        as.data.frame(x@anno)
    } else {
        NULL
    }
})

# Remove NULL entries
annotation_summary <- annotation_summary[!sapply(annotation_summary, is.null)]

# ============================================================
# Save annotations
# ============================================================
cat("\n[", as.character(Sys.time()), "] Saving annotation results...\n")

# Save complete annotation list as RDS
saveRDS(peakAnnoList, file = "annotation/peakAnnoList.rds")
cat("  Saved: annotation/peakAnnoList.rds\n")

# Write summary CSV
annotation_counts <- lapply(annotation_summary, function(df) {
    table(df$annotation)
})
annotation_df <- as.data.frame(do.call(rbind, annotation_counts))
write.csv(annotation_df, "annotation/annotation_summary.csv")
cat("  Saved: annotation/annotation_summary.csv\n")

# Save individual annotation files
for (name in names(annotation_summary)) {
    output_file <- paste0("annotation/", name, "_annotated.csv")
    write.csv(annotation_summary[[name]], output_file, row.names = FALSE)
    cat("  Saved:", output_file, "\n")
}

# ============================================================
# Generate plots
# ============================================================
cat("\n[", as.character(Sys.time()), "] Generating plots...\n")

# Plot 1: Annotation distribution (bar plot)
cat("  Creating annotation bar plot...\n")
pdf("annotation/annotation_barplot.pdf", width = 10, height = 6)
print(plotAnnoBar(peakAnnoList))
dev.off()

png("annotation/annotation_barplot.png", width = 1000, height = 600, res = 100)
print(plotAnnoBar(peakAnnoList))
dev.off()

# Plot 2: Distance to TSS
cat("  Creating distance to TSS plots...\n")
for (name in names(peakAnnoList)) {
    if (!is.null(peakAnnoList[[name]])) {
        pdf(paste0("annotation/", name, "_distance_to_TSS.pdf"), width = 8, height = 6)
        print(plotDistToTSS(peakAnnoList[[name]], title = paste(name, "- Distance to TSS")))
        dev.off()
    }
}

# Plot 3: Annotation pie charts
cat("  Creating annotation pie charts...\n")
for (name in names(peakAnnoList)) {
    if (!is.null(peakAnnoList[[name]])) {
        pdf(paste0("annotation/", name, "_annotation_pie.pdf"), width = 8, height = 8)
        print(plotAnnoPie(peakAnnoList[[name]]))
        dev.off()
    }
}

# ============================================================
# Gene lists
# ============================================================
cat("\n[", as.character(Sys.time()), "] Extracting gene lists...\n")

for (name in names(annotation_summary)) {
    df <- annotation_summary[[name]]

    # Extract unique gene IDs
    genes <- unique(df$geneId)
    genes <- genes[!is.na(genes)]

    # Save gene list
    gene_file <- paste0("annotation/", name, "_genes.txt")
    writeLines(genes, gene_file)
    cat("  Saved:", gene_file, "(", length(genes), "genes )\n")
}

# ============================================================
# Summary report
# ============================================================
cat("\n========================================\n")
cat("Annotation Summary\n")
cat("========================================\n\n")

for (name in names(annotation_summary)) {
    df <- annotation_summary[[name]]
    cat(name, ":\n")
    cat("  Total peaks:", nrow(df), "\n")
    cat("  Unique genes:", length(unique(df$geneId[!is.na(df$geneId)])), "\n")

    # Annotation distribution (top 5)
    anno_counts <- sort(table(df$annotation), decreasing = TRUE)
    cat("  Top annotations:\n")
    for (i in 1:min(5, length(anno_counts))) {
        cat(sprintf("    - %-30s: %d (%.1f%%)\n",
                   names(anno_counts)[i],
                   anno_counts[i],
                   anno_counts[i] * 100 / nrow(df)))
    }
    cat("\n")
}

cat("========================================\n")
cat("Peak Annotation Complete!\n")
cat("========================================\n")
cat("\nOutput files:\n")
cat("  - annotation/peakAnnoList.rds\n")
cat("  - annotation/annotation_summary.csv\n")
cat("  - annotation/*_annotated.csv\n")
cat("  - annotation/*_genes.txt\n")
cat("  - annotation/*.pdf (plots)\n")
cat("========================================\n")
