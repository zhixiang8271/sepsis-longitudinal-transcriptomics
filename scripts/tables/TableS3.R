analysis_file <- file.path("outputs_publication", "analysis.rds")

if (!file.exists(analysis_file)) {
    stop("Required analysis output was not found. Run analysis.R or run.R first.")
}

obj <- readRDS(analysis_file)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

results <- obj$master_results

required_cols <- c(
    "section", "cohort", "role", "interval", "analysis", "specification",
    "effect_measure", "n", "n_observations", "k", "estimate", "SE",
    "CI_low", "CI_high", "p_value", "p_method"
)

missing_cols <- setdiff(required_cols, colnames(results))
if (length(missing_cols) > 0L) {
    stop(
        "analysis_master_results is missing column(s): ",
        paste(missing_cols, collapse = ", ")
    )
}

if (nrow(results) != 26L) {
    stop("Expected 26 rows in analysis_master_results.")
}

fmt_p <- function(p) {
    if (is.na(p)) {
        return("\u2014")
    }

    if (p < 0.001) {
        return("<0.001")
    }

    sprintf("%.3f", p)
}

display_interval <- function(x) {
    out <- x
    out[x == "D1-D2"] <- "D1 \u2192 D2"
    out[x == "D1-D3"] <- "D1 \u2192 D3"
    out[x == "D1-D5"] <- "D1 \u2192 D5"
    out
}

display_rows <- lapply(seq_len(nrow(results)), function(i) {
    row <- results[i, , drop = FALSE]
    is_lmm <- row$effect_measure == "Beta per day"

    if (is_lmm) {
        timepoint <- row$analysis
        n_display <- sprintf(
            "%d patients / %d observations",
            row$n,
            row$n_observations
        )
        effect <- sprintf(
            "\u03b2/day = %.5f (SE %.5f)",
            row$estimate,
            row$SE
        )
    } else {
        timepoint <- display_interval(row$interval)
        n_display <- as.character(row$n)
        effect <- sprintf(
            "\u03c1 = %.3f (95%% CI %.3f to %.3f)",
            row$estimate,
            row$CI_low,
            row$CI_high
        )
    }

    data.frame(
        Section = row$section,
        Cohort = row$cohort,
        `Timepoint / interval` = timepoint,
        Specification = row$specification,
        N = n_display,
        `Effect estimate` = effect,
        P = fmt_p(row$p_value),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
})

table_s3 <- do.call(rbind, display_rows)

write.csv(
    table_s3,
    file.path("tables", "TableS3.csv"),
    row.names = FALSE,
    na = ""
)

write.csv(
    results,
    file.path("tables", "TableS3_numeric.csv"),
    row.names = FALSE,
    na = ""
)

print(table_s3, row.names = FALSE)
