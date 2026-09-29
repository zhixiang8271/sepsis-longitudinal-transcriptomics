get_run_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)

    if (length(file_arg) > 0L) {
        return(sub("^--file=", "", file_arg[1]))
    }

    frame_files <- vapply(
        sys.frames(),
        function(x) {
            if (is.null(x$ofile)) NA_character_ else as.character(x$ofile)
        },
        character(1)
    )

    frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]

    if (length(frame_files) > 0L) {
        return(tail(frame_files, 1))
    }

    "run.R"
}

PROJECT_DIR <- dirname(
    normalizePath(
        get_run_path(),
        winslash = "/",
        mustWork = TRUE
    )
)

setwd(PROJECT_DIR)

required_packages <- c(
    "GEOquery",
    "Biobase",
    "singscore",
    "lme4",
    "lmerTest",
    "emmeans",
    "MCPcounter",
    "illuminaHumanv4.db",
    "AnnotationDbi",
    "dplyr",
    "stringr",
    "ggplot2",
    "patchwork"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
    stop(
        "Missing required R packages: ",
        paste(missing_packages, collapse = ", "),
        "\nInstall the missing packages before running the analysis."
    )
}

run_script <- function(path) {
    full_path <- file.path(PROJECT_DIR, path)

    if (!file.exists(full_path)) {
        stop("Required script was not found: ", path)
    }

    message("Running ", path)
    sys.source(
        full_path,
        envir = new.env(parent = globalenv())
    )
}

run_script("analysis.R")

table_scripts <- c(
    "scripts/tables/TableS1.R",
    "scripts/tables/TableS2.R",
    "scripts/tables/TableS3.R"
)

figure_scripts <- c(
    "scripts/figures/Figure2.R",
    "scripts/figures/Figure3.R",
    "scripts/figures/Figure4.R",
    "scripts/figures/FigureS1.R",
    "scripts/figures/FigureS2.R",
    "scripts/figures/FigureS3.R",
    "scripts/figures/FigureS4.R"
)

invisible(lapply(table_scripts, run_script))
invisible(lapply(figure_scripts, run_script))

message("Complete. Tables are in tables/ and figures are in figures/.")
