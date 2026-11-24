# ==============================================================================
# Project: Differential Gene Expression & Pathway Enrichment in Ythdc1 KO mESCs
# Dataset: NCBI GEO GSE74397 (Roundtree et al., eLife 2017)
# ==============================================================================

library(readxl)
library(ggplot2)
library(ggrepel)
library(gprofiler2)

dir.create("results", showWarnings = FALSE)

# 1. Load & Clean CuffDiff Expression Data
raw_df <- read_excel("data/GSE74397_Ythdc1KO_mESC_RNA-Seq_CuffDiff.xlsx")

clean_df <- subset(raw_df, status == "OK")
clean_df$log2FC <- suppressWarnings(as.numeric(clean_df[["log2(fold_change)"]]))
clean_df$q_val  <- suppressWarnings(as.numeric(clean_df[["q_value"]]))
clean_df <- subset(clean_df, !is.na(log2FC) & !is.na(q_val) & is.finite(log2FC) & q_val > 0)

# 2. Classify Significant DEGs (|log2FC| > 1 & q-value < 0.05)
clean_df$Regulation <- "Not Significant"
clean_df$Regulation[clean_df$log2FC > 1 & clean_df$q_val < 0.05]  <- "Upregulated"
clean_df$Regulation[clean_df$log2FC < -1 & clean_df$q_val < 0.05] <- "Downregulated"

# 3. Map NCBI RefSeq IDs (NM_...) to Official Mouse Gene Symbols (MGI)
clean_df$refseq_clean <- sub("\\..*", "", clean_df$gene_id)
converted <- gconvert(
  query      = unique(clean_df$refseq_clean),
  organism   = "mmusculus",
  target     = "MGI",
  mthreshold = 1,
  filter_na  = TRUE
)
clean_df$Gene_Symbol <- converted$name[match(clean_df$refseq_clean, converted$input)]

sig_genes <- subset(clean_df, Regulation != "Not Significant")
sig_genes <- sig_genes[order(sig_genes$q_val), ]
write.csv(sig_genes, "data/Significant_DEGs_Ythdc1_KO.csv", row.names = FALSE)

# 4. Generate Volcano Plot (Top 10 Annotated Up & Down Genes Labeled)
mapped_only <- subset(clean_df, !is.na(Gene_Symbol) & Gene_Symbol != "None")
mapped_only <- mapped_only[!duplicated(mapped_only$Gene_Symbol), ]

up_sub   <- subset(mapped_only, Regulation == "Upregulated")
top_up   <- head(up_sub[order(up_sub$q_val, -up_sub$log2FC), ], 10)

down_sub <- subset(mapped_only, Regulation == "Downregulated")
top_down <- head(down_sub[order(down_sub$q_val, down_sub$log2FC), ], 10)

top_genes <- rbind(top_up, top_down)

volcano_plot <- ggplot(clean_df, aes(x = log2FC, y = -log10(q_val), color = Regulation)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Downregulated" = "#2b8cbe",
                                "Not Significant" = "grey70",
                                "Upregulated" = "#e34a33")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = Gene_Symbol),
    size         = 3.5,
    fontface     = "bold",
    color        = "black",
    max.overlaps = Inf,
    box.padding  = 0.45
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot") +
  labs(
    title    = "Differential Gene Expression: Ythdc1 KO vs WT (mESCs)",
    subtitle = "Dataset: GSE74397 | Top 10 Annotated DEGs (|log2FC| > 1, q < 0.05)",
    x        = "Log2 Fold Change (KO / WT)",
    y        = "-Log10 Adjusted P-value (q-value)"
  )

ggsave("results/Ythdc1_KO_Volcano_Plot.png", plot = volcano_plot, width = 9.5, height = 6.5, dpi = 300)

# 5. Functional Enrichment Analysis (GO Biological Process & KEGG)
sig_refseq   <- unique(clean_df$refseq_clean[clean_df$Regulation != "Not Significant"])
gost_res     <- gost(query = sig_refseq, organism = "mmusculus", sources = c("GO:BP", "KEGG"))
top_pathways <- head(gost_res$result[order(gost_res$result$p_value), ], 15)

pathway_plot <- ggplot(top_pathways, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value), fill = source)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot") +
  labs(
    title = "Enriched Biological Pathways in Ythdc1 KO (mESCs)",
    x     = NULL,
    y     = "-Log10 Adjusted P-value"
  )

ggsave("results/Ythdc1_KO_Pathway_Enrichment.png", plot = pathway_plot, width = 10, height = 6, dpi = 300)
