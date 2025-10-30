#!/usr/bin/env Rscript
# ============================================================
# Step 7: Functional Enrichment Analysis
# - GO (Gene Ontology) enrichment analysis
# - KEGG pathway enrichment analysis
# - Generate enrichment plots
# ============================================================

cat("========================================\n")
cat("Starting enrichment analysis pipeline\n")
cat("========================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Rn.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
})

# Set working directory to project root
project_dir <- "projects/rat-midbrain"
setwd(project_dir)

# Configuration
anno_dir <- "06_annotation"
output_dir <- "07_enrichment"
results_dir <- "07_enrichment/results"
figures_dir <- "07_enrichment/figures"

# Parameters
pvalue_cutoff <- 0.05
qvalue_cutoff <- 0.2
ontology <- "BP"  # Biological Process

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("Configuration:\n")
cat("  Annotation directory:", anno_dir, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  P-value cutoff:", pvalue_cutoff, "\n")
cat("  Q-value cutoff:", qvalue_cutoff, "\n")
cat("  Ontology:", ontology, "\n\n")

# ============================================================
# Function: Run GO enrichment
# ============================================================
run_go_enrichment <- function(gene_list, sample_name, ont = "BP") {
    cat("\n==========================================\n")
    cat("Running GO enrichment for:", sample_name, "\n")
    cat("Ontology:", ont, "\n")
    cat("==========================================\n")

    if (length(gene_list) == 0) {
        cat("Warning: Empty gene list for", sample_name, "\n")
        return(NULL)
    }

    cat("Input genes:", length(gene_list), "\n")

    tryCatch({
        ego <- enrichGO(
            gene = gene_list,
            OrgDb = org.Rn.eg.db,
            keyType = "ENTREZID",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalue_cutoff,
            qvalueCutoff = qvalue_cutoff,
            readable = TRUE
        )

        if (!is.null(ego) && nrow(ego@result) > 0) {
            cat("  Enriched GO terms:", nrow(ego@result), "\n")
            return(ego)
        } else {
            cat("  No significant GO terms found\n")
            return(NULL)
        }

    }, error = function(e) {
        cat("Error in GO enrichment for", sample_name, ":", conditionMessage(e), "\n")
        return(NULL)
    })
}

# ============================================================
# Function: Run KEGG enrichment
# ============================================================
run_kegg_enrichment <- function(gene_list, sample_name) {
    cat("\n==========================================\n")
    cat("Running KEGG enrichment for:", sample_name, "\n")
    cat("==========================================\n")

    if (length(gene_list) == 0) {
        cat("Warning: Empty gene list for", sample_name, "\n")
        return(NULL)
    }

    cat("Input genes:", length(gene_list), "\n")

    tryCatch({
        kegg <- enrichKEGG(
            gene = gene_list,
            organism = "rno",  # Rattus norvegicus
            pvalueCutoff = pvalue_cutoff,
            qvalueCutoff = qvalue_cutoff
        )

        if (!is.null(kegg) && nrow(kegg@result) > 0) {
            cat("  Enriched KEGG pathways:", nrow(kegg@result), "\n")
            return(kegg)
        } else {
            cat("  No significant KEGG pathways found\n")
            return(NULL)
        }

    }, error = function(e) {
        cat("Error in KEGG enrichment for", sample_name, ":", conditionMessage(e), "\n")
        return(NULL)
    })
}

# ============================================================
# Load gene lists from annotation
# ============================================================
cat("\n[", as.character(Sys.time()), "] Loading gene lists...\n")

gene_files <- list.files(anno_dir, pattern = "_genes\\.txt$", full.names = TRUE)

if (length(gene_files) == 0) {
    cat("Error: No gene files found in", anno_dir, "\n")
    cat("Please run 06_annotation.R first\n")
    quit(status = 1)
}

cat("Found", length(gene_files), "gene list(s)\n\n")

# Load all gene lists
gene_lists <- list()
for (gene_file in gene_files) {
    sample_name <- gsub("_genes\\.txt$", "", basename(gene_file))
    genes <- readLines(gene_file)
    gene_lists[[sample_name]] <- genes
    cat("  Loaded", sample_name, ":", length(genes), "genes\n")
}

# ============================================================
# Run GO enrichment for each sample
# ============================================================
cat("\n[", as.character(Sys.time()), "] Running GO enrichment analyses...\n")

go_results <- list()

for (sample_name in names(gene_lists)) {
    genes <- gene_lists[[sample_name]]

    # Run GO enrichment for BP (Biological Process)
    ego_bp <- run_go_enrichment(genes, sample_name, ont = "BP")
    if (!is.null(ego_bp)) {
        go_results[[paste0(sample_name, "_BP")]] <- ego_bp

        # Save results
        result_file <- file.path(results_dir, paste0(sample_name, "_GO_BP.csv"))
        write.csv(ego_bp@result, result_file, row.names = FALSE)
        cat("  Saved:", result_file, "\n")
    }

    # Run GO enrichment for MF (Molecular Function)
    ego_mf <- run_go_enrichment(genes, sample_name, ont = "MF")
    if (!is.null(ego_mf)) {
        go_results[[paste0(sample_name, "_MF")]] <- ego_mf

        result_file <- file.path(results_dir, paste0(sample_name, "_GO_MF.csv"))
        write.csv(ego_mf@result, result_file, row.names = FALSE)
        cat("  Saved:", result_file, "\n")
    }

    # Run GO enrichment for CC (Cellular Component)
    ego_cc <- run_go_enrichment(genes, sample_name, ont = "CC")
    if (!is.null(ego_cc)) {
        go_results[[paste0(sample_name, "_CC")]] <- ego_cc

        result_file <- file.path(results_dir, paste0(sample_name, "_GO_CC.csv"))
        write.csv(ego_cc@result, result_file, row.names = FALSE)
        cat("  Saved:", result_file, "\n")
    }
}

# ============================================================
# Run KEGG enrichment for each sample
# ============================================================
cat("\n[", as.character(Sys.time()), "] Running KEGG enrichment analyses...\n")

kegg_results <- list()

for (sample_name in names(gene_lists)) {
    genes <- gene_lists[[sample_name]]

    kegg <- run_kegg_enrichment(genes, sample_name)
    if (!is.null(kegg)) {
        kegg_results[[sample_name]] <- kegg

        # Save results
        result_file <- file.path(results_dir, paste0(sample_name, "_KEGG.csv"))
        write.csv(kegg@result, result_file, row.names = FALSE)
        cat("  Saved:", result_file, "\n")
    }
}

# ============================================================
# Generate enrichment plots
# ============================================================
cat("\n[", as.character(Sys.time()), "] Generating enrichment plots...\n")

# GO plots
for (result_name in names(go_results)) {
    ego <- go_results[[result_name]]

    cat("  Creating plots for", result_name, "\n")

    # Dot plot
    pdf(file.path(figures_dir, paste0(result_name, "_dotplot.pdf")), width = 10, height = 8)
    print(dotplot(ego, showCategory = 20, title = paste(result_name, "- Top 20 GO Terms")))
    dev.off()

    png(file.path(figures_dir, paste0(result_name, "_dotplot.png")), width = 1000, height = 800, res = 100)
    print(dotplot(ego, showCategory = 20, title = paste(result_name, "- Top 20 GO Terms")))
    dev.off()

    # Bar plot
    pdf(file.path(figures_dir, paste0(result_name, "_barplot.pdf")), width = 10, height = 8)
    print(barplot(ego, showCategory = 20, title = paste(result_name, "- Top 20 GO Terms")))
    dev.off()

    png(file.path(figures_dir, paste0(result_name, "_barplot.png")), width = 1000, height = 800, res = 100)
    print(barplot(ego, showCategory = 20, title = paste(result_name, "- Top 20 GO Terms")))
    dev.off()

    # Network plot (top 30)
    if (nrow(ego@result) >= 5) {
        tryCatch({
            pdf(file.path(figures_dir, paste0(result_name, "_network.pdf")), width = 12, height = 10)
            print(emapplot(pairwise_termsim(ego), showCategory = 30))
            dev.off()

            png(file.path(figures_dir, paste0(result_name, "_network.png")), width = 1200, height = 1000, res = 100)
            print(emapplot(pairwise_termsim(ego), showCategory = 30))
            dev.off()
        }, error = function(e) {
            cat("    Warning: Could not create network plot:", conditionMessage(e), "\n")
        })
    }
}

# KEGG plots
for (sample_name in names(kegg_results)) {
    kegg <- kegg_results[[sample_name]]

    cat("  Creating KEGG plots for", sample_name, "\n")

    # Dot plot
    pdf(file.path(figures_dir, paste0(sample_name, "_KEGG_dotplot.pdf")), width = 10, height = 8)
    print(dotplot(kegg, showCategory = 20, title = paste(sample_name, "- Top 20 KEGG Pathways")))
    dev.off()

    png(file.path(figures_dir, paste0(sample_name, "_KEGG_dotplot.png")), width = 1000, height = 800, res = 100)
    print(dotplot(kegg, showCategory = 20, title = paste(sample_name, "- Top 20 KEGG Pathways")))
    dev.off()

    # Bar plot
    pdf(file.path(figures_dir, paste0(sample_name, "_KEGG_barplot.pdf")), width = 10, height = 8)
    print(barplot(kegg, showCategory = 20, title = paste(sample_name, "- Top 20 KEGG Pathways")))
    dev.off()

    png(file.path(figures_dir, paste0(sample_name, "_KEGG_barplot.png")), width = 1000, height = 800, res = 100)
    print(barplot(kegg, showCategory = 20, title = paste(sample_name, "- Top 20 KEGG Pathways")))
    dev.off()
}

# ============================================================
# Generate summary report
# ============================================================
cat("\n[", as.character(Sys.time()), "] Generating summary report...\n")

summary_file <- file.path(output_dir, "enrichment_summary.txt")
sink(summary_file)

cat("========================================\n")
cat("Enrichment Analysis Summary\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("========================================\n\n")

cat("Parameters:\n")
cat("  P-value cutoff:", pvalue_cutoff, "\n")
cat("  Q-value cutoff:", qvalue_cutoff, "\n\n")

cat("GO Enrichment Results:\n")
cat("---\n")
for (result_name in names(go_results)) {
    ego <- go_results[[result_name]]
    cat(sprintf("%-30s: %d significant terms\n", result_name, nrow(ego@result)))
}

cat("\nKEGG Enrichment Results:\n")
cat("---\n")
for (sample_name in names(kegg_results)) {
    kegg <- kegg_results[[sample_name]]
    cat(sprintf("%-30s: %d significant pathways\n", sample_name, nrow(kegg@result)))
}

cat("\n========================================\n")

sink()

# Print summary to console
cat("\n")
cat(readLines(summary_file), sep = "\n")

# ============================================================
# Summary
# ============================================================
cat("\n==========================================\n")
cat("Enrichment Analysis Complete!\n")
cat("==========================================\n")
cat("Results:\n")
cat("  - Enrichment results:", results_dir, "\n")
cat("  - Figures:", figures_dir, "\n")
cat("  - Summary:", summary_file, "\n")
cat("\n")
cat("Next step: Run motif analysis (08_motif_analysis.sh)\n")
cat("==========================================\n")
