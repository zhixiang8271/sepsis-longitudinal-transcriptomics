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


scatter_d1 <- bind_rows(
    obj$emtab_d1 %>%
        transmute(
            cohort = "E-MTAB-5273",
            HLA_score = Score_HLA,
            Tcell_score = Score_Tcell
        ),
    obj$gse_d1 %>%
        transmute(
            cohort = "GSE236713",
            HLA_score = Score_HLA,
            Tcell_score = Score_Tcell
        )
)

fixed_timepoint <- obj$fixed_timepoint_raw

COL_PRIMARY <- "#2C7FB8"
COL_VALIDATION <- "#E46C2C"

cohort_colours <- c(
    "E-MTAB-5273" = COL_PRIMARY,
    "GSE236713" = COL_VALIDATION
)

theme_pub <- theme_classic(base_size = 10, base_family = "sans") +
    theme(
        axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.line = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(
            size = 11,
            face = "bold",
            hjust = 0,
            margin = margin(b = 6)
        ),
        legend.position = "none",
        plot.margin = margin(6, 8, 6, 8)
    )

make_d1_scatter <- function(cohort_name, plot_colour) {
    dat <- scatter_d1 %>%
        filter(cohort == cohort_name)

    stat <- fixed_timepoint %>%
        filter(
            cohort == cohort_name,
            timepoint == "D1"
        )

    stopifnot(nrow(stat) == 1L)

    stat_label <- paste0(
        "n = ", stat$n,
        "\nSpearman \u03c1 = ", sprintf("%.2f", stat$rho),
        "\n95% CI ", sprintf("%.2f", stat$CI_low),
        " to ", sprintf("%.2f", stat$CI_high)
    )

    ggplot(
        dat,
        aes(x = HLA_score, y = Tcell_score)
    ) +
        geom_point(
            size = 1.75,
            alpha = 0.63,
            colour = plot_colour
        ) +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linewidth = 0.60,
            colour = plot_colour
        ) +
        annotate(
            "label",
            x = -Inf,
            y = Inf,
            label = stat_label,
            hjust = -0.03,
            vjust = 1.03,
            size = 2.75,
            lineheight = 1.0,
            label.size = 0,
            fill = "white",
            alpha = 0.90
        ) +
        scale_x_continuous(
            expand = expansion(mult = c(0.08, 0.08))
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0.08, 0.10))
        ) +
        labs(
            title = cohort_name,
            x = "HLA-II programme score",
            y = "T-cell signature score"
        ) +
        theme_pub
}

pA <- make_d1_scatter("E-MTAB-5273", COL_PRIMARY)
pB <- make_d1_scatter("GSE236713", COL_VALIDATION)

forest_data <- fixed_timepoint %>%
    mutate(
        cohort = factor(
            cohort,
            levels = c("E-MTAB-5273", "GSE236713")
        ),
        y_pos = case_when(
            cohort == "E-MTAB-5273" & timepoint == "D1" ~ 6,
            cohort == "E-MTAB-5273" & timepoint == "D3" ~ 5,
            cohort == "E-MTAB-5273" & timepoint == "D5" ~ 4,
            cohort == "GSE236713" & timepoint == "D1" ~ 3,
            cohort == "GSE236713" & timepoint == "D2" ~ 2,
            cohort == "GSE236713" & timepoint == "D5" ~ 1
        )
    )

label_data <- forest_data %>%
    arrange(y_pos) %>%
    transmute(
        y_pos,
        label = paste0(timepoint, "  (n=", n, ")")
    )

y_labels <- setNames(
    label_data$label,
    as.character(label_data$y_pos)
)

pC <- ggplot(
    forest_data,
    aes(x = rho, y = y_pos, colour = cohort)
) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4,
        colour = "grey65"
    ) +
    geom_hline(
        yintercept = 3.5,
        linewidth = 0.4,
        colour = "grey82"
    ) +
    geom_errorbarh(
        aes(xmin = CI_low, xmax = CI_high),
        height = 0,
        linewidth = 0.75
    ) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = cohort_colours) +
    scale_y_continuous(
        breaks = 1:6,
        labels = y_labels,
        limits = c(0.6, 6.75),
        expand = c(0, 0)
    ) +
    scale_x_continuous(
        limits = c(0, 0.9),
        breaks = seq(0, 0.8, by = 0.2),
        expand = expansion(mult = c(0.01, 0.02))
    ) +
    annotate(
        "text",
        x = 0.02,
        y = 6.55,
        label = "E-MTAB-5273",
        hjust = 0,
        size = 3.2,
        fontface = "bold",
        colour = COL_PRIMARY
    ) +
    annotate(
        "text",
        x = 0.02,
        y = 3.32,
        label = "GSE236713",
        hjust = 0,
        size = 3.2,
        fontface = "bold",
        colour = COL_VALIDATION
    ) +
    labs(
        title = "Cross-sectional associations across timepoints",
        x = "Spearman \u03c1 (95% CI)",
        y = NULL
    ) +
    theme_pub +
    theme(
        axis.text.y = element_text(size = 9),
        plot.margin = margin(6, 10, 6, 8)
    )

Figure2 <- (
    (pA | pB) /
        patchwork::free(pC, side = "l", type = "space")
) +
    plot_layout(heights = c(1.08, 0.78)) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(size = 12, face = "bold")
    )

print(Figure2)

ggsave(
    file.path("figures", "Figure2.png"),
    Figure2,
    width = 180,
    height = 145,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "Figure2.tiff"),
    Figure2,
    width = 180,
    height = 145,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "Figure2.pdf"),
        Figure2,
        width = 180,
        height = 145,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
