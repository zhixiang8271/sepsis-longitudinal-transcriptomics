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


delta_data <- obj$paired_change_data %>%
    transmute(
        cohort,
        interval,
        patient_id,
        delta_HLA = dHLA,
        delta_Tcell = dTcell
    )

delta_stats <- obj$raw_delta

COL_PRIMARY <- "#2C7FB8"
COL_VALIDATION <- "#E46C2C"

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
            margin = margin(b = 1)
        ),
        plot.subtitle = element_text(
            size = 9.5,
            colour = "black",
            margin = margin(b = 5)
        ),
        legend.position = "none",
        plot.margin = margin(7, 8, 7, 8)
    )

get_limits <- function(dat, cohort_name) {
    d <- dat %>%
        filter(cohort == cohort_name)

    xr <- range(d$delta_HLA, na.rm = TRUE)
    yr <- range(d$delta_Tcell, na.rm = TRUE)

    xpad <- diff(xr) * 0.07
    ypad <- diff(yr) * 0.07

    list(
        x = c(xr[1] - xpad, xr[2] + xpad),
        y = c(yr[1] - ypad, yr[2] + ypad)
    )
}

limits_emtab <- get_limits(delta_data, "E-MTAB-5273")
limits_gse <- get_limits(delta_data, "GSE236713")

make_delta_plot <- function(
    cohort_name,
    interval_name,
    plot_colour,
    axis_limits
) {
    dat <- delta_data %>%
        filter(
            cohort == cohort_name,
            interval == interval_name
        )

    stat <- delta_stats %>%
        filter(
            cohort == cohort_name,
            interval == interval_name
        )

    stopifnot(nrow(stat) == 1L)

    interval_display <- gsub("-", " \u2192 ", interval_name)

    stat_label <- paste0(
        "n = ", stat$n,
        "\nSpearman \u03c1 = ", sprintf("%.2f", stat$rho),
        "\n95% CI ", sprintf("%.2f", stat$CI_low),
        " to ", sprintf("%.2f", stat$CI_high)
    )

    ggplot(
        dat,
        aes(x = delta_HLA, y = delta_Tcell)
    ) +
        geom_hline(
            yintercept = 0,
            linetype = "dashed",
            linewidth = 0.4,
            colour = "grey65"
        ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.4,
            colour = "grey65"
        ) +
        geom_point(
            size = 1.8,
            alpha = 0.65,
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
            size = 2.7,
            lineheight = 1.0,
            label.size = 0,
            fill = "white",
            alpha = 0.90
        ) +
        coord_cartesian(
            xlim = axis_limits$x,
            ylim = axis_limits$y
        ) +
        labs(
            title = interval_display,
            subtitle = cohort_name,
            x = "\u0394 HLA-II programme score",
            y = "\u0394 T-cell signature score"
        ) +
        theme_pub
}

pA <- make_delta_plot(
    "E-MTAB-5273",
    "D1-D3",
    COL_PRIMARY,
    limits_emtab
)

pB <- make_delta_plot(
    "E-MTAB-5273",
    "D1-D5",
    COL_PRIMARY,
    limits_emtab
)

pC <- make_delta_plot(
    "GSE236713",
    "D1-D2",
    COL_VALIDATION,
    limits_gse
)

pD <- make_delta_plot(
    "GSE236713",
    "D1-D5",
    COL_VALIDATION,
    limits_gse
)

Figure3 <- (
    (pA | pB) /
        (pC | pD)
) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(size = 12, face = "bold")
    )

print(Figure3)

ggsave(
    file.path("figures", "Figure3.png"),
    Figure3,
    width = 180,
    height = 150,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "Figure3.tiff"),
    Figure3,
    width = 180,
    height = 150,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "Figure3.pdf"),
        Figure3,
        width = 180,
        height = 150,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
