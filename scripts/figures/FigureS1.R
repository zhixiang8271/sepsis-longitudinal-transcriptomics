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


COL_EMTAB <- "#2C7FB8"
COL_GSE <- "#E46C2C"
COL_TEXT <- "grey15"

availability <- obj$fixed_timepoint_raw %>%
    transmute(
        cohort,
        timepoint,
        n,
        x = case_when(
            timepoint == "D1" ~ 1,
            timepoint %in% c("D2", "D3") ~ 2,
            timepoint == "D5" ~ 3
        )
    )

paired <- obj$raw_delta %>%
    group_by(cohort) %>%
    arrange(interval, .by_group = TRUE) %>%
    mutate(
        start_x = 1,
        end_x = if_else(
            interval %in% c("D1-D2", "D1-D3"),
            2,
            3
        ),
        y = c(2.25, 1.35),
        label = paste0(
            gsub("-", " \u2192 ", interval),
            "    paired n = ",
            n
        )
    ) %>%
    ungroup()

lmm_counts <- obj$lmm_summary %>%
    distinct(cohort, n_patients, n_observations)

theme_pub <- theme_void(base_family = "sans") +
    theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.margin = margin(8, 12, 8, 12)
    )

make_sampling_plot <- function(cohort_name, cohort_colour, panel_tag) {
    dat_avail <- availability %>%
        filter(cohort == cohort_name)

    dat_pair <- paired %>%
        filter(cohort == cohort_name)

    counts <- lmm_counts %>%
        filter(cohort == cohort_name)

    stopifnot(nrow(counts) == 1L)

    ggplot() +
        annotate(
            "text",
            x = 0.65,
            y = 4.75,
            label = "Available samples",
            hjust = 0,
            fontface = "bold",
            size = 3.4
        ) +
        geom_point(
            data = dat_avail,
            aes(x = x, y = 4.00),
            colour = cohort_colour,
            size = 3.8
        ) +
        geom_text(
            data = dat_avail,
            aes(
                x = x,
                y = 3.53,
                label = paste0(timepoint, "\n(n = ", n, ")")
            ),
            size = 3.15,
            lineheight = 1.05
        ) +
        annotate(
            "text",
            x = 0.65,
            y = 2.85,
            label = "Paired longitudinal analyses",
            hjust = 0,
            fontface = "bold",
            size = 3.4
        ) +
        geom_segment(
            data = dat_pair,
            aes(
                x = start_x,
                xend = end_x,
                y = y,
                yend = y
            ),
            colour = cohort_colour,
            linewidth = 0.85
        ) +
        geom_point(
            data = dat_pair,
            aes(x = start_x, y = y),
            colour = cohort_colour,
            size = 2.7
        ) +
        geom_point(
            data = dat_pair,
            aes(x = end_x, y = y),
            colour = cohort_colour,
            size = 2.7
        ) +
        geom_text(
            data = dat_pair,
            aes(
                x = (start_x + end_x) / 2,
                y = y + 0.27,
                label = label
            ),
            size = 2.95
        ) +
        annotate(
            "segment",
            x = 0.65,
            xend = 3.35,
            y = 0.72,
            yend = 0.72,
            colour = "grey82",
            linewidth = 0.55
        ) +
        annotate(
            "text",
            x = 2.00,
            y = 0.42,
            label = paste0(
                "LMM contributors\n",
                counts$n_patients,
                " patients  \u00b7  ",
                counts$n_observations,
                " observations"
            ),
            hjust = 0.5,
            vjust = 0.5,
            fontface = "bold",
            size = 2.85,
            lineheight = 1.15,
            colour = COL_TEXT
        ) +
        coord_cartesian(
            xlim = c(0.55, 3.45),
            ylim = c(0.05, 5.05),
            clip = "off"
        ) +
        labs(
            title = cohort_name,
            tag = panel_tag
        ) +
        theme_pub +
        theme(
            plot.title = element_text(
                colour = cohort_colour,
                size = 12,
                face = "bold"
            ),
            plot.tag.position = c(0, 1)
        )
}

pA <- make_sampling_plot("E-MTAB-5273", COL_EMTAB, "A")
pB <- make_sampling_plot("GSE236713", COL_GSE, "B")

FigureS1 <- pA | pB

print(FigureS1)

ggsave(
    file.path("figures", "FigureS1.png"),
    FigureS1,
    width = 180,
    height = 88,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "FigureS1.tiff"),
    FigureS1,
    width = 180,
    height = 88,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "FigureS1.pdf"),
        FigureS1,
        width = 180,
        height = 88,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
