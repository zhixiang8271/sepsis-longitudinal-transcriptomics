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


primary <- obj$adjusted_primary %>%
    transmute(
        cohort,
        interval,
        specification = "Immune8 (clean markers; primary)",
        n,
        k,
        rho,
        CI_low,
        CI_high
    )

original <- obj$adjusted_original %>%
    transmute(
        cohort,
        interval,
        specification = "Immune8 (original markers)",
        n,
        k,
        rho,
        CI_low,
        CI_high
    )

s4_data <- bind_rows(primary, original) %>%
    mutate(
        specification = factor(
            specification,
            levels = c(
                "Immune8 (clean markers; primary)",
                "Immune8 (original markers)"
            )
        ),
        interval_label = case_when(
            interval == "D1" ~ "D1",
            interval == "D1-D2" ~ "D1 \u2192 D2",
            interval == "D1-D3" ~ "D1 \u2192 D3",
            interval == "D1-D5" ~ "D1 \u2192 D5"
        )
    )

stopifnot(
    nrow(primary) == 6L,
    nrow(original) == 6L,
    !anyNA(s4_data$interval_label),
    !anyNA(s4_data$rho),
    !anyNA(s4_data$CI_low),
    !anyNA(s4_data$CI_high)
)

COL_PRIMARY <- "#2C7FB8"
COL_ORIGINAL <- "#8C6BB1"

spec_colours <- c(
    "Immune8 (clean markers; primary)" = COL_PRIMARY,
    "Immune8 (original markers)" = COL_ORIGINAL
)

spec_shapes <- c(
    "Immune8 (clean markers; primary)" = 15,
    "Immune8 (original markers)" = 17
)

theme_s4 <- theme_classic(base_size = 10, base_family = "sans") +
    theme(
        axis.title = element_text(size = 9.7, colour = "black"),
        axis.text = element_text(size = 8.9, colour = "black"),
        axis.line = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(
            size = 11,
            face = "bold",
            hjust = 0,
            margin = margin(b = 6)
        ),
        plot.tag = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8.7),
        legend.position = "bottom",
        plot.margin = margin(8, 10, 8, 10)
    )

make_forest <- function(
    cohort_name,
    interval_order,
    cohort_colour,
    panel_tag
) {
    interval_order_display <- case_when(
        interval_order == "D1" ~ "D1",
        interval_order == "D1-D2" ~ "D1 \u2192 D2",
        interval_order == "D1-D3" ~ "D1 \u2192 D3",
        interval_order == "D1-D5" ~ "D1 \u2192 D5"
    )

    dat <- s4_data %>%
        filter(cohort == cohort_name) %>%
        mutate(
            interval_label = factor(
                interval_label,
                levels = rev(interval_order_display)
            )
        )

    dodge <- position_dodge(width = 0.38)

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
            linewidth = 0.70,
            position = dodge
        ) +
        geom_point(
            size = 2.7,
            position = dodge
        ) +
        scale_colour_manual(
            values = spec_colours,
            breaks = c(
                "Immune8 (clean markers; primary)",
                "Immune8 (original markers)"
            ),
            drop = FALSE
        ) +
        scale_shape_manual(
            values = spec_shapes,
            breaks = c(
                "Immune8 (clean markers; primary)",
                "Immune8 (original markers)"
            ),
            drop = FALSE
        ) +
        scale_x_continuous(
            limits = c(-0.8, 0.8),
            breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
            labels = c("-0.8", "-0.4", "0", "0.4", "0.8"),
            expand = c(0, 0)
        ) +
        labs(
            tag = panel_tag,
            title = cohort_name,
            x = "Partial Spearman \u03c1 (95% CI)",
            y = NULL,
            colour = NULL,
            shape = NULL
        ) +
        theme_s4 +
        theme(
            plot.title = element_text(colour = cohort_colour)
        )
}

pA <- make_forest(
    "E-MTAB-5273",
    c("D1", "D1-D3", "D1-D5"),
    "#2C7FB8",
    "A"
)

pB <- make_forest(
    "GSE236713",
    c("D1", "D1-D2", "D1-D5"),
    "#E46C2C",
    "B"
)

FigureS4 <- (
    pA | pB
) +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center"
    )

print(FigureS4)

ggsave(
    file.path("figures", "FigureS4.png"),
    FigureS4,
    width = 190,
    height = 100,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "FigureS4.tiff"),
    FigureS4,
    width = 190,
    height = 100,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "FigureS4.pdf"),
        FigureS4,
        width = 190,
        height = 100,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
