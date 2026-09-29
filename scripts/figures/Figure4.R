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


raw_d1 <- obj$raw_d1 %>%
    transmute(
        cohort,
        interval = timepoint,
        specification = "Raw",
        rho,
        CI_low,
        CI_high
    )

raw_delta <- obj$raw_delta %>%
    transmute(
        cohort,
        interval,
        specification = "Raw",
        rho,
        CI_low,
        CI_high
    )

adjusted <- obj$adjusted_primary %>%
    transmute(
        cohort,
        interval,
        specification = "Immune8-adjusted",
        rho,
        CI_low,
        CI_high
    )

fig4 <- bind_rows(raw_d1, raw_delta, adjusted) %>%
    mutate(
        interval_label = case_when(
            interval == "D1" ~ "D1",
            interval == "D1-D2" ~ "D1 \u2192 D2",
            interval == "D1-D3" ~ "D1 \u2192 D3",
            interval == "D1-D5" ~ "D1 \u2192 D5"
        ),
        specification = factor(
            specification,
            levels = c("Raw", "Immune8-adjusted")
        )
    )

COL_RAW <- "grey25"
COL_ADJ <- "#2C7FB8"

spec_cols <- c(
    "Raw" = COL_RAW,
    "Immune8-adjusted" = COL_ADJ
)

spec_shapes <- c(
    "Raw" = 16,
    "Immune8-adjusted" = 15
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
            margin = margin(b = 8)
        ),
        plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(8, 10, 7, 8)
    )

make_forest <- function(cohort_name, interval_order, panel_tag) {
    interval_order_display <- case_when(
        interval_order == "D1" ~ "D1",
        interval_order == "D1-D2" ~ "D1 \u2192 D2",
        interval_order == "D1-D3" ~ "D1 \u2192 D3",
        interval_order == "D1-D5" ~ "D1 \u2192 D5"
    )

    dat <- fig4 %>%
        filter(cohort == cohort_name) %>%
        mutate(
            interval_label = factor(
                interval_label,
                levels = rev(interval_order_display)
            )
        )

    dodge <- position_dodge(width = 0.48)

    ggplot(
        dat,
        aes(
            x = rho,
            y = interval_label,
            colour = specification,
            shape = specification
        )
    ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.4,
            colour = "grey65"
        ) +
        geom_errorbarh(
            aes(xmin = CI_low, xmax = CI_high),
            height = 0,
            linewidth = 0.75,
            position = dodge
        ) +
        geom_point(
            size = 2.8,
            position = dodge
        ) +
        scale_colour_manual(values = spec_cols) +
        scale_shape_manual(values = spec_shapes) +
        scale_x_continuous(
            limits = c(-0.65, 0.90),
            breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9),
            labels = c("-0.6", "-0.3", "0", "0.3", "0.6", "0.9"),
            expand = c(0, 0)
        ) +
        labs(
            title = cohort_name,
            tag = panel_tag,
            x = "Association estimate (\u03c1, 95% CI)",
            y = NULL
        ) +
        theme_pub +
        theme(
            plot.tag.position = c(0, 1),
            plot.tag = element_text(
                size = 12,
                face = "bold",
                hjust = 0,
                vjust = 1
            )
        )
}

pA <- make_forest(
    "E-MTAB-5273",
    c("D1", "D1-D3", "D1-D5"),
    "A"
)

pB <- make_forest(
    "GSE236713",
    c("D1", "D1-D2", "D1-D5"),
    "B"
)

legend_df <- data.frame(
    x = c(1.8, 4.9),
    y = c(1.0, 1.0),
    label = c("Raw", "Immune8-adjusted"),
    colour = c(COL_RAW, COL_ADJ),
    shape = c(16, 15)
)

legend_plot <- ggplot() +
    geom_segment(
        data = legend_df,
        aes(
            x = x - 0.25,
            xend = x + 0.25,
            y = y,
            yend = y
        ),
        colour = legend_df$colour,
        linewidth = 0.75
    ) +
    geom_point(
        data = legend_df,
        aes(x = x, y = y),
        colour = legend_df$colour,
        shape = legend_df$shape,
        size = 2.8
    ) +
    geom_text(
        data = legend_df,
        aes(
            x = x + 0.38,
            y = y,
            label = label
        ),
        hjust = 0,
        vjust = 0.5,
        size = 3.0,
        family = "sans"
    ) +
    scale_x_continuous(
        limits = c(0, 8.0),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        limits = c(0.65, 1.35),
        expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
        plot.margin = margin(0, 8, 0, 8)
    )

Figure4 <- (pA | pB) /
    legend_plot +
    plot_layout(heights = c(1, 0.10))

print(Figure4)

ggsave(
    file.path("figures", "Figure4.png"),
    Figure4,
    width = 190,
    height = 112,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "Figure4.tiff"),
    Figure4,
    width = 190,
    height = 112,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "Figure4.pdf"),
        Figure4,
        width = 190,
        height = 112,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
