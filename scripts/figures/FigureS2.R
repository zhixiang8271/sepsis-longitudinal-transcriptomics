suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

analysis_file <- file.path("outputs_publication", "analysis.rds")

if (!file.exists(analysis_file)) {
    stop(
        "Required analysis output was not found. ",
        "Run analysis.R or run.R from the repository root first."
    )
}

obj <- readRDS(analysis_file)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)


hla_genes <- obj$hla_19
tcell_genes <- obj$tcell_14

gse_expr <- obj$gse_d1_expression
emtab_expr <- obj$emtab_d1_expression

stopifnot(
    all(hla_genes %in% rownames(gse_expr)),
    all(hla_genes %in% rownames(emtab_expr)),
    all(tcell_genes %in% rownames(gse_expr)),
    all(tcell_genes %in% rownames(emtab_expr))
)

gse_scores <- setNames(
    obj$gse_d1$Score_HLA,
    obj$gse_d1$geo_accession
)

emtab_scores <- setNames(
    obj$emtab_d1$Score_HLA,
    obj$emtab_d1$sample_id
)

gse_tcell_scores <- setNames(
    obj$gse_d1$Score_Tcell,
    obj$gse_d1$geo_accession
)

emtab_tcell_scores <- setNames(
    obj$emtab_d1$Score_Tcell,
    obj$emtab_d1$sample_id
)

run_pca <- function(expr, score_lookup, cohort) {
    x <- t(expr)

    keep <- apply(
        x,
        2,
        function(z) {
            s <- sd(z, na.rm = TRUE)
            is.finite(s) && s > 0
        }
    )

    x <- x[, keep, drop = FALSE]

    pca <- prcomp(
        x,
        center = TRUE,
        scale. = TRUE
    )

    pc1_scores <- pca$x[, "PC1"]
    pc1_loadings <- pca$rotation[, "PC1"]

    r <- suppressWarnings(
        cor(
            pc1_scores,
            score_lookup[rownames(x)],
            method = "spearman",
            use = "complete.obs"
        )
    )

    if (is.finite(r) && r < 0) {
        pc1_scores <- -pc1_scores
        pc1_loadings <- -pc1_loadings
    }

    data.frame(
        cohort = cohort,
        gene = names(pc1_loadings),
        PC1_loading = as.numeric(pc1_loadings),
        stringsAsFactors = FALSE
    )
}

make_loading_comparison <- function(
    gse_expression,
    emtab_expression,
    gse_score,
    emtab_score
) {
    gse_loading <- run_pca(
        gse_expression,
        gse_score,
        "GSE236713"
    )

    emtab_loading <- run_pca(
        emtab_expression,
        emtab_score,
        "E-MTAB-5273"
    )

    names(gse_loading)[3] <- "GSE_PC1_loading"
    names(emtab_loading)[3] <- "EMTAB_PC1_loading"

    out <- merge(
        gse_loading[, c("gene", "GSE_PC1_loading")],
        emtab_loading[, c("gene", "EMTAB_PC1_loading")],
        by = "gene"
    )

    out$same_direction <- sign(out$GSE_PC1_loading) ==
        sign(out$EMTAB_PC1_loading)

    out
}

hla_gse <- gse_expr[hla_genes, , drop = FALSE]
hla_emtab <- emtab_expr[hla_genes, , drop = FALSE]
tcell_gse <- gse_expr[tcell_genes, , drop = FALSE]
tcell_emtab <- emtab_expr[tcell_genes, , drop = FALSE]

hla_cor_gse <- cor(
    t(hla_gse),
    method = "spearman",
    use = "pairwise.complete.obs"
)

hla_cor_emtab <- cor(
    t(hla_emtab),
    method = "spearman",
    use = "pairwise.complete.obs"
)

tcell_cor_gse <- cor(
    t(tcell_gse),
    method = "spearman",
    use = "pairwise.complete.obs"
)

tcell_cor_emtab <- cor(
    t(tcell_emtab),
    method = "spearman",
    use = "pairwise.complete.obs"
)

hla_loadings <- make_loading_comparison(
    hla_gse,
    hla_emtab,
    gse_scores,
    emtab_scores
)

tcell_loadings <- make_loading_comparison(
    tcell_gse,
    tcell_emtab,
    gse_tcell_scores,
    emtab_tcell_scores
)

extract_cor <- function(mat, cohort, construct) {
    mat <- as.matrix(mat)

    data.frame(
        rho = mat[upper.tri(mat)],
        cohort = cohort,
        construct = construct,
        stringsAsFactors = FALSE
    )
}

cor_data <- bind_rows(
    extract_cor(hla_cor_emtab, "E-MTAB-5273", "HLA-II programme"),
    extract_cor(hla_cor_gse, "GSE236713", "HLA-II programme"),
    extract_cor(tcell_cor_emtab, "E-MTAB-5273", "T-cell signature"),
    extract_cor(tcell_cor_gse, "GSE236713", "T-cell signature")
)

hla_loading_r <- cor(
    hla_loadings$GSE_PC1_loading,
    hla_loadings$EMTAB_PC1_loading,
    method = "spearman",
    use = "complete.obs"
)

tcell_loading_r <- cor(
    tcell_loadings$GSE_PC1_loading,
    tcell_loadings$EMTAB_PC1_loading,
    method = "spearman",
    use = "complete.obs"
)

hla_same_dir <- sum(
    hla_loadings$same_direction,
    na.rm = TRUE
)

tcell_same_dir <- sum(
    tcell_loadings$same_direction,
    na.rm = TRUE
)

make_square_limits <- function(x, y, padding = 0.08) {
    vals <- c(x, y)
    vals <- vals[is.finite(vals)]
    lim <- range(vals)
    span <- diff(lim)

    c(
        lim[1] - span * padding,
        lim[2] + span * padding
    )
}

hla_lim <- make_square_limits(
    hla_loadings$GSE_PC1_loading,
    hla_loadings$EMTAB_PC1_loading
)

tcell_lim <- make_square_limits(
    tcell_loadings$GSE_PC1_loading,
    tcell_loadings$EMTAB_PC1_loading
)

COL_EMTAB <- "#2C7FB8"
COL_GSE <- "#E46C2C"
COL_HLA <- "#4C78A8"
COL_TCELL <- "#9B6FB6"

theme_top <- theme_classic(base_size = 10, base_family = "sans") +
    theme(
        axis.title = element_text(size = 9.5, colour = "black"),
        axis.text = element_text(size = 8.8, colour = "black"),
        axis.line = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0,
            margin = margin(b = 5)
        ),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(-0.08, 1.04),
        legend.position = "none",
        plot.margin = margin(8, 12, 8, 12)
    )

theme_bottom <- theme_classic(base_size = 10, base_family = "sans") +
    theme(
        axis.title = element_text(size = 9.3, colour = "black"),
        axis.text = element_text(size = 8.7, colour = "black"),
        axis.line = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0,
            margin = margin(b = 2)
        ),
        plot.subtitle = element_text(
            size = 8.2,
            colour = "grey25",
            hjust = 0,
            lineheight = 1.12,
            margin = margin(b = 7)
        ),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(-0.08, 1.04),
        aspect.ratio = 1,
        legend.position = "none",
        plot.margin = margin(9, 16, 12, 16)
    )

pA_data <- cor_data %>%
    filter(construct == "HLA-II programme") %>%
    mutate(
        cohort = factor(
            cohort,
            levels = c("E-MTAB-5273", "GSE236713")
        )
    )

pA <- ggplot(
    pA_data,
    aes(x = cohort, y = rho, fill = cohort)
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4,
        colour = "grey70"
    ) +
    geom_violin(
        width = 0.75,
        alpha = 0.18,
        linewidth = 0.55
    ) +
    geom_boxplot(
        width = 0.20,
        outlier.shape = NA,
        linewidth = 0.55,
        alpha = 0.80
    ) +
    scale_fill_manual(
        values = c(
            "E-MTAB-5273" = COL_EMTAB,
            "GSE236713" = COL_GSE
        )
    ) +
    scale_y_continuous(
        limits = c(-1, 1),
        breaks = seq(-1, 1, 0.5)
    ) +
    labs(
        tag = "A",
        title = "HLA-II programme",
        x = NULL,
        y = "Pairwise gene\u2013gene Spearman \u03c1"
    ) +
    theme_top

pB_data <- cor_data %>%
    filter(construct == "T-cell signature") %>%
    mutate(
        cohort = factor(
            cohort,
            levels = c("E-MTAB-5273", "GSE236713")
        )
    )

pB <- ggplot(
    pB_data,
    aes(x = cohort, y = rho, fill = cohort)
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4,
        colour = "grey70"
    ) +
    geom_violin(
        width = 0.75,
        alpha = 0.18,
        linewidth = 0.55
    ) +
    geom_boxplot(
        width = 0.20,
        outlier.shape = NA,
        linewidth = 0.55,
        alpha = 0.80
    ) +
    scale_fill_manual(
        values = c(
            "E-MTAB-5273" = COL_EMTAB,
            "GSE236713" = COL_GSE
        )
    ) +
    scale_y_continuous(
        limits = c(-1, 1),
        breaks = seq(-1, 1, 0.5)
    ) +
    labs(
        tag = "B",
        title = "T-cell signature",
        x = NULL,
        y = "Pairwise gene\u2013gene Spearman \u03c1"
    ) +
    theme_top

pC <- ggplot(
    hla_loadings,
    aes(
        x = GSE_PC1_loading,
        y = EMTAB_PC1_loading
    )
) +
    geom_hline(
        yintercept = 0,
        linewidth = 0.4,
        colour = "grey82"
    ) +
    geom_vline(
        xintercept = 0,
        linewidth = 0.4,
        colour = "grey82"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.45,
        colour = "grey65"
    ) +
    geom_point(
        size = 2.8,
        alpha = 0.85,
        colour = COL_HLA
    ) +
    scale_x_continuous(
        limits = hla_lim,
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        limits = hla_lim,
        expand = c(0, 0)
    ) +
    labs(
        tag = "C",
        title = "HLA-II programme",
        subtitle = paste0(
            "Spearman \u03c1 = ",
            sprintf("%.2f", hla_loading_r),
            "\nConcordant loadings ",
            hla_same_dir,
            "/",
            nrow(hla_loadings)
        ),
        x = "GSE236713 PC1 loading",
        y = "E-MTAB-5273 PC1 loading"
    ) +
    theme_bottom

pD <- ggplot(
    tcell_loadings,
    aes(
        x = GSE_PC1_loading,
        y = EMTAB_PC1_loading
    )
) +
    geom_hline(
        yintercept = 0,
        linewidth = 0.4,
        colour = "grey82"
    ) +
    geom_vline(
        xintercept = 0,
        linewidth = 0.4,
        colour = "grey82"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.45,
        colour = "grey65"
    ) +
    geom_point(
        size = 2.8,
        alpha = 0.85,
        colour = COL_TCELL
    ) +
    scale_x_continuous(
        limits = tcell_lim,
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        limits = tcell_lim,
        expand = c(0, 0)
    ) +
    labs(
        tag = "D",
        title = "T-cell signature",
        subtitle = paste0(
            "Spearman \u03c1 = ",
            sprintf("%.2f", tcell_loading_r),
            "\nConcordant loadings ",
            tcell_same_dir,
            "/",
            nrow(tcell_loadings)
        ),
        x = "GSE236713 PC1 loading",
        y = "E-MTAB-5273 PC1 loading"
    ) +
    theme_bottom

pc1_header <- wrap_elements(
    full = grid::textGrob(
        "PC1 loading agreement across cohorts",
        gp = grid::gpar(
            fontfamily = "sans",
            fontsize = 11,
            fontface = "bold"
        )
    )
)

top_row <- pA | pB

bottom_panels <- (
    pC | plot_spacer() | pD
) +
    plot_layout(widths = c(1, 0.12, 1))

bottom_row <- (
    pc1_header /
        bottom_panels
) +
    plot_layout(heights = c(0.08, 1))

FigureS2 <- (
    top_row /
        bottom_row
) +
    plot_layout(heights = c(0.88, 1.18))

print(FigureS2)

ggsave(
    file.path("figures", "FigureS2.png"),
    FigureS2,
    width = 200,
    height = 165,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "FigureS2.tiff"),
    FigureS2,
    width = 200,
    height = 165,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "FigureS2.pdf"),
        FigureS2,
        width = 200,
        height = 165,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
