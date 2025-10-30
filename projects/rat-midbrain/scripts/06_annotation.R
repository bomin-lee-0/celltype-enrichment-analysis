#!/usr/bin/env Rscript
# ============================================================
# Step 6: Peak Annotation with ChIPseeker
# - Annotate peaks to genomic features
# - Generate annotation plots
# - Export annotated peak tables
# ============================================================

cat("========================================\n")
cat("Starting peak annotation pipeline\n")
cat("========================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
    library(org.Rn.eg.db)
    library(GenomicRanges)
    library(ggplot2)
    library(dplyr)
})

# Set working directory to project root
project_dir <- "projects/rat-midbrain"
setwd(project_dir)

# Configuration
peak_dir <- "05_peak_processing/unique_peaks"
merged_dir <- "05_peak_processing/merged_peaks"
output_dir <- "06_annotation"
plot_dir <- "06_annotation/plots"

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load TxDb for rat genome
txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene

# TSS region definition
tss_region <- c(-3000, 3000)

cat("Configuration:\n")
cat("  Peak directory:", peak_dir, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  TSS region:", tss_region[1], "to", tss_region[2], "bp\n\n")

# ============================================================
# Function: Annotate peaks
# ============================================================
annotate_peaks <- function(peak_file, sample_name, txdb, tss_region) {
    cat("\n==========================================\n")
    cat("Annotating:", sample_name, "\n")
    cat("==========================================\n")

    tryCatch({
        # Read and annotate peaks
        peak_anno <- annotatePeak(
            peak_file,
            TxDb = txdb,
            tssRegion = tss_region,
            verbose = FALSE
        )

        # Print summary
        cat("\nAnnotation summary:\n")
        print(peak_anno)

        # Convert to data frame
        peak_anno_df <- as.data.frame(peak_anno)

        # Add sample name
        peak_anno_df$sample <- sample_name

        cat("  Total peaks annotated:", nrow(peak_anno_df), "\n")

        return(list(
            annotation = peak_anno,
            dataframe = peak_anno_df
        ))

    }, error = function(e) {
        cat("Error annotating", sample_name, ":", conditionMessage(e), "\n")
        return(NULL)
    })
}

# ============================================================
# Annotate unique peaks for each cell type
# ============================================================
cat("\n[", as.character(Sys.time()), "] Annotating unique peaks...\n")

# Find all unique peak files
unique_peak_files <- list.files(peak_dir, pattern = "_unique\\.bed$", full.names = TRUE)

if (length(unique_peak_files) == 0) {
    cat("Warning: No unique peak files found in", peak_dir, "\n")
    cat("Checking merged peaks instead...\n")
    unique_peak_files <- list.files(merged_dir, pattern = "_merged\\.bed$", full.names = TRUE)
}

cat("Found", length(unique_peak_files), "peak file(s) to annotate\n\n")

# Store annotations
all_annotations <- list()
all_dataframes <- list()

for (peak_file in unique_peak_files) {
    # Extract sample name
    sample_name <- gsub("_unique\\.bed$|_merged\\.bed$", "", basename(peak_file))

    # Annotate
    result <- annotate_peaks(peak_file, sample_name, txdb, tss_region)

    if (!is.null(result)) {
        all_annotations[[sample_name]] <- result$annotation
        all_dataframes[[sample_name]] <- result$dataframe

        # Save individual annotation file
        output_file <- file.path(output_dir, paste0(sample_name, "_annotated.csv"))
        write.csv(result$dataframe, output_file, row.names = FALSE)
        cat("  Saved:", output_file, "\n")
    }
}

# ============================================================
# Generate annotation plots
# ============================================================
cat("\n[", as.character(Sys.time()), "] Generating annotation plots...\n")

if (length(all_annotations) > 0) {

    # Plot 1: Annotation pie chart for each sample
    for (sample_name in names(all_annotations)) {
        cat("  Creating pie chart for", sample_name, "\n")

        pdf(file.path(plot_dir, paste0(sample_name, "_annotation_pie.pdf")), width = 8, height = 8)
        print(plotAnnoPie(all_annotations[[sample_name]]))
        dev.off()

        png(file.path(plot_dir, paste0(sample_name, "_annotation_pie.png")), width = 800, height = 800, res = 100)
        print(plotAnnoPie(all_annotations[[sample_name]]))
        dev.off()
    }

    # Plot 2: Annotation bar plot for each sample
    for (sample_name in names(all_annotations)) {
        cat("  Creating bar plot for", sample_name, "\n")

        pdf(file.path(plot_dir, paste0(sample_name, "_annotation_bar.pdf")), width = 10, height = 6)
        print(plotAnnoBar(all_annotations[[sample_name]]))
        dev.off()

        png(file.path(plot_dir, paste0(sample_name, "_annotation_bar.png")), width = 1000, height = 600, res = 100)
        print(plotAnnoBar(all_annotations[[sample_name]]))
        dev.off()
    }

    # Plot 3: Distance to TSS distribution
    for (sample_name in names(all_annotations)) {
        cat("  Creating TSS distance plot for", sample_name, "\n")

        pdf(file.path(plot_dir, paste0(sample_name, "_tss_distance.pdf")), width = 10, height = 6)
        print(plotDistToTSS(all_annotations[[sample_name]], title = paste(sample_name, "- Distance to TSS")))
        dev.off()

        png(file.path(plot_dir, paste0(sample_name, "_tss_distance.png")), width = 1000, height = 600, res = 100)
        print(plotDistToTSS(all_annotations[[sample_name]], title = paste(sample_name, "- Distance to TSS")))
        dev.off()
    }

    # Plot 4: Compare annotations across samples (if multiple)
    if (length(all_annotations) > 1) {
        cat("  Creating comparison plots\n")

        # Create named list for comparison
        anno_list <- all_annotations

        pdf(file.path(plot_dir, "annotation_comparison.pdf"), width = 12, height = 8)
        print(plotAnnoBar(anno_list))
        dev.off()

        png(file.path(plot_dir, "annotation_comparison.png"), width = 1200, height = 800, res = 100)
        print(plotAnnoBar(anno_list))
        dev.off()
    }
}

# ============================================================
# Generate annotation summary
# ============================================================
cat("\n[", as.character(Sys.time()), "] Generating annotation summary...\n")

if (length(all_dataframes) > 0) {
    # Combine all annotations
    combined_df <- bind_rows(all_dataframes)

    # Summary statistics
    summary_stats <- combined_df %>%
        group_by(sample, annotation) %>%
        summarise(count = n(), .groups = "drop") %>%
        arrange(sample, desc(count))

    # Save summary
    summary_file <- file.path(output_dir, "annotation_summary.csv")
    write.csv(summary_stats, summary_file, row.names = FALSE)
    cat("  Saved summary:", summary_file, "\n")

    # Print summary
    cat("\nAnnotation Summary by Sample:\n")
    print(summary_stats)

    # Gene-level summary
    gene_summary <- combined_df %>%
        filter(!is.na(geneId)) %>%
        group_by(sample) %>%
        summarise(
            total_peaks = n(),
            unique_genes = n_distinct(geneId),
            .groups = "drop"
        )

    cat("\nGene-level Summary:\n")
    print(gene_summary)

    # Save gene list for each sample
    for (sample_name in names(all_dataframes)) {
        genes <- all_dataframes[[sample_name]] %>%
            filter(!is.na(geneId)) %>%
            pull(geneId) %>%
            unique()

        gene_file <- file.path(output_dir, paste0(sample_name, "_genes.txt"))
        writeLines(genes, gene_file)
        cat("  Saved gene list:", gene_file, "(", length(genes), "genes )\n")
    }
}

# ============================================================
# Summary
# ============================================================
cat("\n==========================================\n")
cat("Peak Annotation Complete!\n")
cat("==========================================\n")
cat("Results:\n")
cat("  - Annotated peaks:", output_dir, "\n")
cat("  - Plots:", plot_dir, "\n")
cat("  - Summary:", file.path(output_dir, "annotation_summary.csv"), "\n")
cat("\n")
cat("Next step: Run enrichment analysis (07_enrichment.R)\n")
cat("==========================================\n")
