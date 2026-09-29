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


trajectory <- obj$lmm_trajectories
lmm_summary <- obj$lmm_summary

raw_temporal <- bind_rows(
    obj$emtab_long %>%
        transmute(
            cohort = "E-MTAB-5273",
            day = time_num,
            module = "HLA-II programme",
            score = Score_HLA
        ),
    obj$emtab_long %>%
        transmute(
            cohort = "E-MTAB-5273",
            day = time_num,
            module = "T-cell signature",
            score = Score_Tcell
        ),
    obj$gse_long %>%
        transmute(
            cohort = "GSE236713",
            day = time_num,
            module = "HLA-II programme",
            score = Score_HLA
        ),
    obj$gse_long %>%
        transmute(
            cohort = "GSE236713",
            day = time_num,
            module = "T-cell signature",
            score = Score_Tcell
        )
)

standardise_construct <- function(x) {
    case_when(
        grepl("HLA", x, ignore.case = TRUE) ~ "HLA-II programme",
        grepl("T.?cell|Tcell|dysfunction", x, ignore.case = TRUE) ~ "T-cell signature",
        TRUE ~ x
    )
}

trajectory <- trajectory %>%
    mutate(module_short = standardise_construct(construct))

lmm_summary <- lmm_summary %>%
    mutate(module_short = standardise_construct(construct))

raw_temporal <- raw_temporal %>%
    mutate(module_short = standardise_construct(module))

COL_EMTAB <- "#2C7FB8"
COL_GSE <- "#E46C2C"

fmt_p <- function(p) {
    if (is.na(p)) {
        return("NA")
    }

    if (p < 0.001) {
        return("< 0.001")
    }

    paste0("= ", sprintf("%.3f", p))
}

theme_s3 <- theme_classic(base_size = 10, base_family = "sans") +
    theme(
        axis.title = element_text(size = 9.7, colour = "black"),
        axis.text = element_text(size = 8.9, colour = "black"),
        axis.line = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(
            size = 10.8,
            face = "bold",
            hjust = 0,
            margin = margin(b = 2)
        ),
        plot.subtitle = element_text(
            size = 8.7,
            colour = "grey20",
            hjust = 0,
            lineheight = 1.10,
            margin = margin(b = 5)
        ),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(-0.07, 1.03),
        legend.position = "none",
        plot.margin = margin(8, 12, 8, 12)
    )

make_lmm_plot <- function(
    cohort_name,
    module_name,
    colour,
    panel_tag
) {
    raw_dat <- raw_temporal %>%
        filter(
            cohort == cohort_name,
            module_short == module_name
        )

    traj_dat <- trajectory %>%
        filter(
            cohort == cohort_name,
            module_short == module_name
        ) %>%
        arrange(day)

    stat_dat <- lmm_summary %>%
        filter(
            cohort == cohort_name,
            module_short == module_name
        )

    if (nrow(stat_dat) != 1L) {
        stop(
            "Expected one LMM result for ",
            cohort_name,
            " / ",
            module_name,
            "."
        )
    }

    subtitle_text <- paste0(
        module_name,
        "\n\u03b2/day = ",
        sprintf("%.4f", stat_dat$slope_per_day),
        "   |   P ",
        fmt_p(stat_dat$p_value)
    )

    ggplot() +
        geom_jitter(
            data = raw_dat,
            aes(x = day, y = score),
            width = 0.07,
            height = 0,
            size = 1.0,
            alpha = 0.12,
            colour = colour
        ) +
        geom_errorbar(
            data = traj_dat,
            aes(
                x = day,
                ymin = CI_low,
                ymax = CI_high
            ),
            width = 0.10,
            linewidth = 0.75,
            colour = colour
        ) +
        geom_line(
            data = traj_dat,
            aes(
                x = day,
                y = estimated_score,
                group = 1
            ),
            linewidth = 0.90,
            colour = colour
        ) +
        geom_point(
            data = traj_dat,
            aes(x = day, y = estimated_score),
            size = 2.9,
            colour = colour
        ) +
        scale_x_continuous(
            breaks = traj_dat$day,
            labels = traj_dat$timepoint,
            expand = expansion(mult = c(0.08, 0.08))
        ) +
        labs(
            tag = panel_tag,
            title = cohort_name,
            subtitle = subtitle_text,
            x = "Sampling timepoint",
            y = "singscore"
        ) +
        theme_s3
}

pA <- make_lmm_plot(
    "E-MTAB-5273",
    "HLA-II programme",
    COL_EMTAB,
    "A"
)

pB <- make_lmm_plot(
    "GSE236713",
    "HLA-II programme",
    COL_GSE,
    "B"
)

pC <- make_lmm_plot(
    "E-MTAB-5273",
    "T-cell signature",
    COL_EMTAB,
    "C"
)

pD <- make_lmm_plot(
    "GSE236713",
    "T-cell signature",
    COL_GSE,
    "D"
)

FigureS3 <- (pA | pB) / (pC | pD)

print(FigureS3)

ggsave(
    file.path("figures", "FigureS3.png"),
    FigureS3,
    width = 190,
    height = 155,
    units = "mm",
    dpi = 300
)

ggsave(
    file.path("figures", "FigureS3.tiff"),
    FigureS3,
    width = 190,
    height = 155,
    units = "mm",
    dpi = 600,
    compression = "lzw"
)

if (capabilities("cairo")) {
    ggsave(
        file.path("figures", "FigureS3.pdf"),
        FigureS3,
        width = 190,
        height = 155,
        units = "mm",
        device = cairo_pdf
    )
} else {
    message("Cairo support is unavailable; PDF output was skipped.")
}
