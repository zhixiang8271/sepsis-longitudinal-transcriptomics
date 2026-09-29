# Longitudinal whole-blood transcriptomic analysis in sepsis
# E-MTAB-5273 (primary cohort) and GSE236713 (validation cohort)
# R analysis script accompanying the manuscript

# Configuration
PROJECT_DIR <- normalizePath(".", winslash = "/", mustWork = FALSE)

EMTAB_LOCAL_DIRS <- unique(Filter(
    nzchar,
    c(
        Sys.getenv("EMTAB_DIR", unset = ""),
        PROJECT_DIR
    )
))

GEO_LOCAL_DIRS <- unique(Filter(
    nzchar,
    c(
        Sys.getenv("GEO_LOCAL_DIR", unset = ""),
        PROJECT_DIR
    )
))

OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs_publication")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Reproducibility settings
N_BOOT <- 5000L
PARTIAL_PERM_N <- 5000L

FIXED_BOOT_SEEDS <- c(
    "E-MTAB-5273_D1" = 20260826L,
    "E-MTAB-5273_D3" = 20260827L,
    "E-MTAB-5273_D5" = 20260828L,
    "GSE236713_D1" = 20260825L,
    "GSE236713_D2" = 20260829L,
    "GSE236713_D5" = 20260830L
)

DELTA_BOOT_SEEDS <- c(
    "E-MTAB-5273_D1-D3" = 20260835L,
    "E-MTAB-5273_D1-D5" = 20260836L,
    "GSE236713_D1-D2" = 20260837L,
    "GSE236713_D1-D5" = 20260838L
)

ADJUSTED_BOOT_SEEDS <- c(
    "E-MTAB-5273_D1" = 20260926L,
    "GSE236713_D1" = 20260928L,
    "E-MTAB-5273_D1-D3" = 20261026L,
    "E-MTAB-5273_D1-D5" = 20261028L,
    "GSE236713_D1-D2" = 20261030L,
    "GSE236713_D1-D5" = 20261032L
)

ADJUSTED_PERM_SEEDS <- c(
    "E-MTAB-5273_D1" = 20260915L,
    "GSE236713_D1" = 20260914L,
    "E-MTAB-5273_D1-D3" = 20261016L,
    "E-MTAB-5273_D1-D5" = 20261017L,
    "GSE236713_D1-D2" = 20261014L,
    "GSE236713_D1-D5" = 20261015L
)

SENSITIVITY_BOOT_SEEDS <- c(
    "E-MTAB-5273_D1" = 202609087L,
    "E-MTAB-5273_D1-D3" = 202609088L,
    "E-MTAB-5273_D1-D5" = 202609089L,
    "GSE236713_D1" = 202609090L,
    "GSE236713_D1-D2" = 202609091L,
    "GSE236713_D1-D5" = 202609092L
)

# Network settings
DOWNLOAD_TIMEOUT_SEC <- 600L
DOWNLOAD_MAX_RETRIES <- 3L
DOWNLOAD_RETRY_WAIT  <- 5L

options(timeout = DOWNLOAD_TIMEOUT_SEC)
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = DOWNLOAD_TIMEOUT_SEC)

# Packages
required_packages <- c(
    "GEOquery",
    "singscore",
    "lme4",
    "lmerTest",
    "emmeans",
    "MCPcounter",
    "illuminaHumanv4.db",
    "AnnotationDbi",
    "dplyr",
    "stringr"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
    stop(
        "Missing required R packages:\n  ",
        paste(missing_packages, collapse = ", "),
        "\nInstall them before running this script. ",
        "Bioconductor packages should be installed with BiocManager."
    )
}

suppressPackageStartupMessages({
    library(GEOquery)
    library(singscore)
    library(lme4)
    library(lmerTest)
    library(emmeans)
    library(MCPcounter)
    library(illuminaHumanv4.db)
    library(AnnotationDbi)
    library(dplyr)
    library(stringr)
})

# Gene sets
# 19-gene HLA-II antigen-presentation transcriptional programme
hla_19 <- c(
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
    "CD74", "CIITA", "RFX5", "RFXAP", "RFXANK",
    "CTSS", "LGMN", "IFI30"
)

# 14-gene T-cell dysfunction-associated transcriptional signature
tcell_14 <- c(
    "HAVCR2", "LAG3", "TIGIT", "CTLA4", "BTLA", "CD160",
    "TOX2", "NR4A1", "NR4A2", "BATF", "IRF4",
    "PRDM1", "EOMES", "ENTPD1"
)

stopifnot(length(hla_19) == 19L, length(tcell_14) == 14L)

# MCPcounter Immune8 markers
MCP_MARKERS_IMMUNE8 <- data.frame(
    `HUGO symbols` = c(
        # T cells (16)
        "CD28","CD3D","CD3G","CD5","CD6","CHRM3-AS2","CTLA4","FLT3LG",
        "ICOS","MAL","MGC40069","PBX4","SIRPG","THEMIS","TNFRSF25","TRAT1",
        # CD8 T cells (1)
        "CD8B",
        # Cytotoxic lymphocytes (7)
        "CD8A","EOMES","FGFBP2","GNLY","KLRC3","KLRC4","KLRD1",
        # B lineage (9)
        "BANK1","CD19","CD22","CD79A","CR2","FCRL2","IGKC","MS4A1","PAX5",
        # NK cells (9)
        "CD160","KIR2DL1","KIR2DL3","KIR2DL4","KIR3DL1","KIR3DS1",
        "NCR1","PTGDR","SH2D1B",
        # Monocytic lineage (7)
        "ADAP2","CSF1R","FPR3","KYNU","PLA2G7","RASSF4","TFEC",
        # Myeloid dendritic cells (6)
        "CD1A","CD1B","CD1E","CLEC10A","CLIC2","WFDC21P",
        # Neutrophils (15)
        "CA4","CEACAM3","CXCR1","CXCR2","CYP4F3","FCGR3B","HAL","KCNJ15",
        "MEGF9","SLC25A37","STEAP4","TECPR2","TLE3","TNFRSF10C","VNN3"
    ),
    `Cell population` = c(
        rep("T cells", 16),
        rep("CD8 T cells", 1),
        rep("Cytotoxic lymphocytes", 7),
        rep("B lineage", 9),
        rep("NK cells", 9),
        rep("Monocytic lineage", 7),
        rep("Myeloid dendritic cells", 6),
        rep("Neutrophils", 15)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

stopifnot(nrow(MCP_MARKERS_IMMUNE8) == 70L)

IMMUNE8 <- c(
    "T cells",
    "CD8 T cells",
    "Cytotoxic lymphocytes",
    "B lineage",
    "NK cells",
    "Monocytic lineage",
    "Myeloid dendritic cells",
    "Neutrophils"
)

overlap_with_constructs <- intersect(
    MCP_MARKERS_IMMUNE8$`HUGO symbols`,
    union(hla_19, tcell_14)
)

MCP_MARKERS_CLEAN <- MCP_MARKERS_IMMUNE8[
    !MCP_MARKERS_IMMUNE8$`HUGO symbols` %in% union(hla_19, tcell_14),
    ,
    drop = FALSE
]

expected_overlap <- sort(c("CTLA4", "EOMES", "CD160"))
if (!identical(sort(overlap_with_constructs), expected_overlap)) {
    stop(
        "Unexpected MCPcounter/study-construct overlap. Observed: ",
        paste(sort(overlap_with_constructs), collapse = ", ")
    )
}

# Helper functions
collapse_probes_by_symbol <- function(expr_probe, symbols) {
    keep <- !is.na(symbols) & symbols != ""
    expr_probe <- expr_probe[keep, , drop = FALSE]
    symbols <- symbols[keep]

    split_idx <- split(seq_along(symbols), symbols)

    out <- vapply(
        split_idx,
        function(ii) {
            if (length(ii) == 1L) {
                as.numeric(expr_probe[ii, , drop = TRUE])
            } else {
                apply(expr_probe[ii, , drop = FALSE], 2, median, na.rm = TRUE)
            }
        },
        numeric(ncol(expr_probe))
    )

    out <- t(out)
    colnames(out) <- colnames(expr_probe)
    out
}

score_constructs <- function(expr_gene, hla_genes = hla_19, tcell_genes = tcell_14) {
    missing_hla <- setdiff(hla_genes, rownames(expr_gene))
    missing_tcell <- setdiff(tcell_genes, rownames(expr_gene))

    if (length(missing_hla) > 0L || length(missing_tcell) > 0L) {
        stop(
            "Gene sets are not fully represented in this matrix.\n",
            "Missing HLA-II genes: ",
            ifelse(length(missing_hla) == 0L, "none", paste(missing_hla, collapse = ", ")),
            "\nMissing T-cell signature genes: ",
            ifelse(length(missing_tcell) == 0L, "none", paste(missing_tcell, collapse = ", "))
        )
    }

    ranked <- singscore::rankGenes(expr_gene)

    hla_score <- singscore::simpleScore(
        ranked,
        upSet = hla_genes,
        knownDirection = TRUE
    )$TotalScore

    tcell_score <- singscore::simpleScore(
        ranked,
        upSet = tcell_genes,
        knownDirection = TRUE
    )$TotalScore

    data.frame(
        Sample = colnames(expr_gene),
        Score_HLA = as.numeric(hla_score),
        Score_Tcell = as.numeric(tcell_score),
        stringsAsFactors = FALSE
    )
}

mcp_estimate_local <- function(expression, markers) {
    features <- markers[
        markers$`HUGO symbols` %in% rownames(expression),
        ,
        drop = FALSE
    ]

    features_split <- split(
        features$`HUGO symbols`,
        features$`Cell population`
    )

    # appendSignatures returns one score per signature and sample;
    # transpose to population x sample orientation used below.
    out <- t(
        MCPcounter::appendSignatures(
            expression,
            features_split
        )
    )

    out
}

spearman_test <- function(x, y) {
    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 4L) {
        return(list(rho = NA_real_, p = NA_real_, n = length(x)))
    }

    ct <- suppressWarnings(
        cor.test(x, y, method = "spearman", exact = FALSE)
    )

    list(
        rho = unname(ct$estimate),
        p = ct$p.value,
        n = length(x)
    )
}

bootstrap_spearman_ci <- function(x, y, B = N_BOOT, seed) {
    dat <- data.frame(x = x, y = y)
    dat <- dat[complete.cases(dat), , drop = FALSE]

    n <- nrow(dat)
    if (n < 4L) {
        return(c(rho = NA_real_, lower = NA_real_, upper = NA_real_))
    }

    point <- unname(cor(dat$x, dat$y, method = "spearman"))

    set.seed(seed)

    boot_rho <- replicate(
        B,
        {
            ii <- sample.int(n, n, replace = TRUE)
            suppressWarnings(
                cor(
                    dat$x[ii],
                    dat$y[ii],
                    method = "spearman",
                    use = "complete.obs"
                )
            )
        }
    )

    ci <- quantile(
        boot_rho,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
    )

    c(rho = point, lower = ci[1], upper = ci[2])
}

partial_spearman_rank <- function(x, y, Z) {
    Z <- as.data.frame(Z, check.names = FALSE)

    dat <- data.frame(
        x = as.numeric(x),
        y = as.numeric(y),
        Z,
        check.names = FALSE
    )

    dat <- dat[complete.cases(dat), , drop = FALSE]

    if (nrow(dat) < 5L) {
        stop("Too few complete observations for partial Spearman correlation.")
    }

    cov_names <- colnames(dat)[-(1:2)]

    # Rank-based partial Spearman adjustment:
    # numeric conversion, removal of zero/non-finite-variance covariates,
    # rank transformation, and QR-based removal of exact linear dependencies.
    if (length(cov_names) > 0L) {
        for (nm in cov_names) {
            dat[[nm]] <- as.numeric(dat[[nm]])
        }

        keep_var <- vapply(
            dat[cov_names],
            function(v) {
                vv <- stats::var(v)
                is.finite(vv) && vv > 0
            },
            logical(1)
        )

        cov_names <- cov_names[keep_var]
    }

    rank_vec <- function(v) {
        rank(v, ties.method = "average", na.last = "keep")
    }

    rx <- rank_vec(dat$x)
    ry <- rank_vec(dat$y)

    if (length(cov_names) > 0L) {
        rz <- as.data.frame(
            lapply(dat[cov_names], rank_vec),
            check.names = FALSE
        )

        X0 <- cbind(`(Intercept)` = 1, as.matrix(rz))
        qx <- qr(X0)
        keep_idx <- sort(qx$pivot[seq_len(qx$rank)])
        X <- X0[, keep_idx, drop = FALSE]
    } else {
        X <- matrix(
            1,
            nrow = nrow(dat),
            ncol = 1L,
            dimnames = list(NULL, "(Intercept)")
        )
    }

    fit_x <- lm.fit(X, rx)
    fit_y <- lm.fit(X, ry)

    rho <- cor(
        fit_x$residuals,
        fit_y$residuals,
        method = "pearson"
    )

    list(
        rho = unname(rho),
        n = nrow(dat),
        k = ncol(X) - 1L,
        X = X,
        residual_x = fit_x$residuals,
        residual_y = fit_y$residuals,
        fitted_y = fit_y$fitted.values
    )
}

freedman_lane_p <- function(x, y, Z, nperm = PARTIAL_PERM_N, seed) {
    obs <- partial_spearman_rank(x, y, Z)

    if (nperm <= 0L) {
        return(NA_real_)
    }

    set.seed(seed)

    perm_rho <- numeric(nperm)

    for (b in seq_len(nperm)) {
        y_perm <- obs$fitted_y + sample(
            obs$residual_y,
            length(obs$residual_y),
            replace = FALSE
        )

        residual_y_perm <- lm.fit(
            obs$X,
            y_perm
        )$residuals

        perm_rho[b] <- cor(
            obs$residual_x,
            residual_y_perm,
            method = "pearson"
        )
    }

    (sum(abs(perm_rho) >= abs(obs$rho), na.rm = TRUE) + 1) /
        (sum(is.finite(perm_rho)) + 1)
}

bootstrap_partial_ci <- function(x, y, Z, B = N_BOOT, seed) {
    Z <- as.data.frame(Z, check.names = FALSE)

    dat <- data.frame(
        x = as.numeric(x),
        y = as.numeric(y),
        Z,
        check.names = FALSE
    )

    dat <- dat[complete.cases(dat), , drop = FALSE]

    point <- partial_spearman_rank(
        dat$x,
        dat$y,
        dat[, -(1:2), drop = FALSE]
    )$rho

    n <- nrow(dat)
    set.seed(seed)

    boot_rho <- replicate(
        B,
        {
            ii <- sample.int(n, n, replace = TRUE)
            d <- dat[ii, , drop = FALSE]

            res <- tryCatch(
                partial_spearman_rank(
                    d$x,
                    d$y,
                    d[, -(1:2), drop = FALSE]
                ),
                error = function(e) NULL
            )

            if (is.null(res)) NA_real_ else res$rho
        }
    )

    ci <- quantile(
        boot_rho,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
    )

    c(rho = point, lower = ci[1], upper = ci[2])
}

get_mcp_covariates <- function(mcp_matrix, sample_ids, populations = IMMUNE8) {
    available <- intersect(populations, rownames(mcp_matrix))

    if (length(available) != length(populations)) {
        stop(
            "Not all Immune8 populations are available. Missing: ",
            paste(setdiff(populations, available), collapse = ", ")
        )
    }

    missing_samples <- setdiff(sample_ids, colnames(mcp_matrix))
    if (length(missing_samples) > 0L) {
        stop(
            "MCPcounter matrix is missing sample(s): ",
            paste(head(missing_samples, 10), collapse = ", ")
        )
    }

    out <- as.data.frame(
        t(mcp_matrix[available, sample_ids, drop = FALSE]),
        check.names = FALSE
    )

    colnames(out) <- make.names(colnames(out), unique = TRUE)
    rownames(out) <- sample_ids

    out
}

make_paired_changes <- function(
    long_df,
    sample_col,
    followup_day,
    cohort,
    interval
) {
    d1 <- long_df[
        long_df$time_num == 1,
        c("patient_id", sample_col, "Score_HLA", "Score_Tcell"),
        drop = FALSE
    ]

    followup <- long_df[
        long_df$time_num == followup_day,
        c("patient_id", sample_col, "Score_HLA", "Score_Tcell"),
        drop = FALSE
    ]

    names(d1) <- c("patient_id", "sample_D1", "HLA_D1", "Tcell_D1")
    names(followup) <- c("patient_id", "sample_FU", "HLA_FU", "Tcell_FU")

    pair <- merge(d1, followup, by = "patient_id")

    data.frame(
        cohort = cohort,
        interval = interval,
        patient_id = pair$patient_id,
        sample_D1 = pair$sample_D1,
        sample_FU = pair$sample_FU,
        HLA_D1 = pair$HLA_D1,
        HLA_FU = pair$HLA_FU,
        Tcell_D1 = pair$Tcell_D1,
        Tcell_FU = pair$Tcell_FU,
        dHLA = pair$HLA_FU - pair$HLA_D1,
        dTcell = pair$Tcell_FU - pair$Tcell_D1,
        stringsAsFactors = FALSE
    )
}

make_delta_mcp <- function(pair_df, mcp_matrix, populations = IMMUNE8) {
    cov_d1 <- get_mcp_covariates(
        mcp_matrix,
        pair_df$sample_D1,
        populations
    )

    cov_fu <- get_mcp_covariates(
        mcp_matrix,
        pair_df$sample_FU,
        populations
    )

    cov_d1 <- cov_d1[pair_df$sample_D1, , drop = FALSE]
    cov_fu <- cov_fu[pair_df$sample_FU, , drop = FALSE]

    out <- cov_fu
    for (j in seq_len(ncol(out))) {
        out[[j]] <- cov_fu[[j]] - cov_d1[[j]]
    }

    rownames(out) <- pair_df$patient_id
    out
}

run_adjusted_analysis <- function(x, y, Z, boot_seed, perm_seed) {
    point <- partial_spearman_rank(x, y, Z)

    ci <- bootstrap_partial_ci(
        x,
        y,
        Z,
        B = N_BOOT,
        seed = boot_seed
    )

    p_perm <- freedman_lane_p(
        x,
        y,
        Z,
        nperm = PARTIAL_PERM_N,
        seed = perm_seed
    )

    data.frame(
        n = point$n,
        k = point$k,
        rho = point$rho,
        CI_low = unname(ci["lower"]),
        CI_high = unname(ci["upper"]),
        p_value = p_perm,
        stringsAsFactors = FALSE
    )
}

run_sensitivity_analysis <- function(x, y, Z, boot_seed) {
    point <- partial_spearman_rank(x, y, Z)

    ci <- bootstrap_partial_ci(
        x,
        y,
        Z,
        B = N_BOOT,
        seed = boot_seed
    )

    data.frame(
        n = point$n,
        k = point$k,
        rho = point$rho,
        CI_low = unname(ci["lower"]),
        CI_high = unname(ci["upper"]),
        p_value = NA_real_,
        stringsAsFactors = FALSE
    )
}

# Data access
first_existing_file <- function(directories, filenames) {
    for (dir_i in directories) {
        if (is.na(dir_i) || !nzchar(dir_i) || !dir.exists(dir_i)) {
            next
        }

        for (name_i in filenames) {
            candidate <- file.path(dir_i, name_i)
            if (file.exists(candidate)) {
                return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
            }
        }
    }

    NA_character_
}

download_with_retry <- function(
    url,
    destfile,
    max_retries = DOWNLOAD_MAX_RETRIES,
    min_bytes = 100L
) {
    dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)

    # Use a temporary file so interrupted downloads do not overwrite the cache.
    tmpfile <- paste0(destfile, ".part")
    if (file.exists(tmpfile)) {
        unlink(tmpfile)
    }

    for (attempt in seq_len(max_retries)) {
        cat(
            sprintf(
                "  Online attempt %d/%d: %s\n",
                attempt,
                max_retries,
                url
            )
        )

        ok <- tryCatch(
            {
                suppressWarnings(
                    utils::download.file(
                        url = url,
                        destfile = tmpfile,
                        mode = "wb",
                        quiet = TRUE,
                        method = "libcurl"
                    )
                )

                file.exists(tmpfile) &&
                    is.finite(file.info(tmpfile)$size) &&
                    file.info(tmpfile)$size >= min_bytes
            },
            error = function(e) {
                cat("    Download error: ", conditionMessage(e), "\n", sep = "")
                FALSE
            },
            warning = function(w) {
                cat("    Download warning: ", conditionMessage(w), "\n", sep = "")
                FALSE
            }
        )

        if (isTRUE(ok)) {
            if (file.exists(destfile)) {
                unlink(destfile)
            }

            renamed <- file.rename(tmpfile, destfile)
            if (!renamed) {
                file.copy(tmpfile, destfile, overwrite = TRUE)
                unlink(tmpfile)
            }

            cat(
                "  Online download succeeded: ",
                normalizePath(destfile, winslash = "/", mustWork = TRUE),
                "\n",
                sep = ""
            )
            return(normalizePath(destfile, winslash = "/", mustWork = TRUE))
        }

        if (file.exists(tmpfile)) {
            unlink(tmpfile)
        }

        if (attempt < max_retries) {
            Sys.sleep(DOWNLOAD_RETRY_WAIT)
        }
    }

    NA_character_
}

obtain_public_file <- function(
    label,
    url,
    local_dirs,
    local_filenames,
    download_file,
    min_bytes = 100L
) {
    cat("\n", label, "\n", sep = "")

    local_file <- first_existing_file(
        local_dirs,
        local_filenames
    )

    if (!is.na(local_file)) {
        file_size <- file.info(local_file)$size

        if (is.finite(file_size) && file_size >= min_bytes) {
            cat("Using local file: ", local_file, "\n", sep = "")
            return(local_file)
        }

        cat(
            "Local file is incomplete or unexpectedly small; ",
            "downloading a fresh copy.\n",
            sep = ""
        )
    } else {
        cat("Local file not found; downloading from the public repository.\n")
    }

    online_file <- download_with_retry(
        url = url,
        destfile = download_file,
        min_bytes = min_bytes
    )

    if (!is.na(online_file)) {
        return(online_file)
    }

    stop(
        label,
        " could not be found locally or downloaded.\n",
        "Place one of the following files in:\n  ",
        PROJECT_DIR,
        "\nExpected filename(s):\n  ",
        paste(local_filenames, collapse = "\n  ")
    )
}

getGEO_with_retry <- function(
    geo_id,
    max_retries = DOWNLOAD_MAX_RETRIES,
    ...
) {
    for (attempt in seq_len(max_retries)) {
        cat(
            sprintf(
                "  GEO online attempt %d/%d: %s\n",
                attempt,
                max_retries,
                geo_id
            )
        )

        obj <- tryCatch(
            GEOquery::getGEO(geo_id, ...),
            error = function(e) {
                cat(
                    "    GEO error: ",
                    conditionMessage(e),
                    "\n",
                    sep = ""
                )
                NULL
            }
        )

        if (!is.null(obj)) {
            cat("  GEO online retrieval succeeded: ", geo_id, "\n", sep = "")
            return(obj)
        }

        if (attempt < max_retries) {
            Sys.sleep(DOWNLOAD_RETRY_WAIT)
        }
    }

    NULL
}

# GSE236713
load_gse236713 <- function() {
    cat("\nGSE236713 Series Matrix\n")

    local_series <- first_existing_file(
        GEO_LOCAL_DIRS,
        c(
            "GSE236713_series_matrix.txt.gz",
            "GSE236713_series_matrix.txt"
        )
    )

    gse <- NULL

    if (!is.na(local_series)) {
        cat("Using local file: ", local_series, "\n", sep = "")
        gse <- tryCatch(
            GEOquery::getGEO(filename = local_series),
            error = function(e) {
                stop(
                    "Local GSE236713 Series Matrix could not be parsed: ",
                    conditionMessage(e)
                )
            }
        )
    } else {
        cat("Local file not found; downloading from NCBI GEO.\n")

        series_file <- download_with_retry(
            url = paste0(
                "https://ftp.ncbi.nlm.nih.gov/geo/series/",
                "GSE236nnn/GSE236713/matrix/",
                "GSE236713_series_matrix.txt.gz"
            ),
            destfile = file.path(
                PROJECT_DIR,
                "GSE236713_series_matrix.txt.gz"
            ),
            min_bytes = 10000L
        )

        if (!is.na(series_file)) {
            gse <- tryCatch(
                GEOquery::getGEO(filename = series_file),
                error = function(e) {
                    stop(
                        "Downloaded GSE236713 Series Matrix could not be parsed: ",
                        conditionMessage(e)
                    )
                }
            )
        } else {
            cat("Direct download failed; trying GEOquery retrieval.\n")

            gse <- getGEO_with_retry(
                "GSE236713",
                GSEMatrix = TRUE,
                getGPL = FALSE,
                destdir = PROJECT_DIR
            )
        }
    }

    if (is.null(gse)) {
        stop(
            "GSE236713 could not be found locally or downloaded.\n",
            "Place GSE236713_series_matrix.txt.gz in:\n  ",
            PROJECT_DIR
        )
    }

    if (is.list(gse) && !inherits(gse, "ExpressionSet")) {
        eset <- gse[[1]]
    } else {
        eset <- gse
    }

    list(
        expr_raw = Biobase::exprs(eset),
        meta = Biobase::pData(eset)
    )
}

standardise_gpl17077_annotation <- function(tab) {
    normalise_colname <- function(x) {
        gsub("[^a-z0-9]+", "", tolower(x))
    }

    nm <- colnames(tab)
    nm_norm <- normalise_colname(nm)

    id_candidates <- c(
        "id",
        "idref",
        "probeid",
        "probeidentifier"
    )

    symbol_candidates <- c(
        "genesymbol",
        "genesymbols",
        "symbol"
    )

    id_idx <- match(id_candidates, nm_norm, nomatch = 0L)
    id_idx <- id_idx[id_idx > 0L]

    symbol_idx <- match(symbol_candidates, nm_norm, nomatch = 0L)
    symbol_idx <- symbol_idx[symbol_idx > 0L]

    if (length(id_idx) == 0L || length(symbol_idx) == 0L) {
        stop(
            "GPL17077 annotation does not contain identifiable probe-ID and gene-symbol columns.\n",
            "Columns found: ",
            paste(nm, collapse = ", ")
        )
    }

    data.frame(
        ID = as.character(tab[[id_idx[1]]]),
        GENE_SYMBOL = as.character(tab[[symbol_idx[1]]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

load_gpl17077_table <- function() {
    cat("\nGPL17077 annotation\n")

    local_annot <- first_existing_file(
        GEO_LOCAL_DIRS,
        c(
            "GPL17077.annot.gz",
            "GPL17077.annot",
            "GPL17077.txt"
        )
    )

    if (is.na(local_annot)) {
        for (dir_i in GEO_LOCAL_DIRS) {
            if (is.na(dir_i) || !nzchar(dir_i) || !dir.exists(dir_i)) {
                next
            }

            candidates <- list.files(
                dir_i,
                pattern = "^GPL17077-[0-9]+\\.txt(\\.gz)?$",
                full.names = TRUE,
                ignore.case = TRUE
            )

            if (length(candidates) > 0L) {
                local_annot <- normalizePath(
                    sort(candidates)[1],
                    winslash = "/",
                    mustWork = TRUE
                )
                break
            }
        }
    }

    if (!is.na(local_annot)) {
        cat("Using local file: ", local_annot, "\n", sep = "")

        if (grepl("\\.gz$", local_annot, ignore.case = TRUE)) {
            con <- gzfile(local_annot, open = "rt")
            on.exit(close(con), add = TRUE)
            tab <- read.delim(
                con,
                comment.char = "#",
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        } else {
            tab <- read.delim(
                local_annot,
                comment.char = "#",
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        }

        return(
            standardise_gpl17077_annotation(tab)
        )
    }

    local_soft <- first_existing_file(
        GEO_LOCAL_DIRS,
        c(
            "GPL17077_family.soft.gz",
            "GPL17077_family.soft"
        )
    )

    if (!is.na(local_soft)) {
        cat("Using local file: ", local_soft, "\n", sep = "")

        gpl_local <- GEOquery::getGEO(
            filename = local_soft
        )

        return(
            standardise_gpl17077_annotation(
                GEOquery::Table(gpl_local)
            )
        )
    }

    cat("Local file not found; downloading the GEO platform SOFT file.\n")

    soft_file <- download_with_retry(
        url = paste0(
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
            "GPL17nnn/GPL17077/soft/GPL17077_family.soft.gz"
        ),
        destfile = file.path(
            PROJECT_DIR,
            "GPL17077_family.soft.gz"
        ),
        min_bytes = 1000L
    )

    if (!is.na(soft_file)) {
        gpl_local <- tryCatch(
            GEOquery::getGEO(filename = soft_file),
            error = function(e) NULL
        )

        if (!is.null(gpl_local)) {
            return(
                standardise_gpl17077_annotation(
                    GEOquery::Table(gpl_local)
                )
            )
        }
    }

    cat("Direct downloads failed; trying GEOquery retrieval.\n")

    gpl <- getGEO_with_retry(
        "GPL17077",
        destdir = PROJECT_DIR
    )

    if (!is.null(gpl)) {
        return(
            standardise_gpl17077_annotation(
                GEOquery::Table(gpl)
            )
        )
    }

    stop(
        "GPL17077 could not be found locally or downloaded.\n",
        "Place a GPL17077 annotation file or GPL17077_family.soft.gz in the same directory as analysis.R:\n  ",
        PROJECT_DIR
    )
}

gse_obj <- load_gse236713()
expr_raw_gse <- gse_obj$expr_raw
meta_all_gse <- gse_obj$meta

meta_all_gse$group <- ifelse(
    grepl("Sepsis", meta_all_gse$title),
    "Sepsis",
    ifelse(
        grepl("^Healthy Control_", meta_all_gse$title),
        "HC",
        ifelse(
            grepl("^OOHCA SIRS_", meta_all_gse$title),
            "OOHCA_SIRS",
            "Unknown"
        )
    )
)

meta_gse <- meta_all_gse[
    meta_all_gse$group == "Sepsis",
    ,
    drop = FALSE
]

meta_gse$timepoint <- stringr::str_extract(
    meta_gse$title,
    "Day: ([^,]+)",
    group = 1
)
meta_gse$timepoint <- trimws(meta_gse$timepoint)
meta_gse$timepoint <- ifelse(
    meta_gse$timepoint == "Discharge from ICU",
    "Discharge",
    paste0("D", meta_gse$timepoint)
)

meta_gse$patient_id <- stringr::str_extract(
    meta_gse$title,
    "Patient (\\d+)",
    group = 1
)

meta_gse$infection_source <- ifelse(
    grepl("^Pulmonary", meta_gse$title),
    "Pulmonary",
    ifelse(
        grepl("^Abdominal", meta_gse$title),
        "Abdominal",
        NA_character_
    )
)

outcome_columns <- names(meta_gse)[
    vapply(
        meta_gse,
        function(x) {
            any(
                grepl(
                    "died/survived:",
                    as.character(x),
                    ignore.case = TRUE
                ),
                na.rm = TRUE
            )
        },
        logical(1)
    )
]

if (length(outcome_columns) != 1L) {
    stop("Could not identify a unique 28-day outcome field in GSE236713 metadata.")
}

outcome_raw <- stringr::str_extract(
    as.character(meta_gse[[outcome_columns[1]]]),
    "died/survived: (.+)",
    group = 1
)
outcome_raw <- tolower(trimws(outcome_raw))

meta_gse$outcome_28d <- ifelse(
    outcome_raw == "survived",
    "Survived",
    ifelse(
        outcome_raw == "died",
        "Died",
        NA_character_
    )
)

# Prefer the unsuffixed sample when a repeated patient-timepoint is deposited.
dup_key <- paste(meta_gse$patient_id, meta_gse$timepoint, sep = "__")
dup_keys <- unique(
    dup_key[duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE)]
)

if (length(dup_keys) > 0L) {
    keep <- rep(TRUE, nrow(meta_gse))

    for (dk in dup_keys) {
        ii <- which(dup_key == dk)
        titles <- meta_gse$title[ii]
        is_repeat <- grepl("Patient [0-9]+-[0-9]+", titles)

        if (sum(!is_repeat) == 1L) {
            keep[ii[is_repeat]] <- FALSE
        } else {
            stop(
                "Unable to adjudicate duplicated GSE236713 patient-timepoint: ",
                dk
            )
        }
    }

    meta_gse <- meta_gse[keep, , drop = FALSE]
}

rownames(meta_gse) <- meta_gse$geo_accession

gpl_table <- load_gpl17077_table()

if (!all(c("ID", "GENE_SYMBOL") %in% colnames(gpl_table))) {
    stop("GPL17077 annotation does not contain expected ID and GENE_SYMBOL columns.")
}

probe_to_gene <- data.frame(
    probe_id = gpl_table$ID,
    symbol = gpl_table$GENE_SYMBOL,
    stringsAsFactors = FALSE
)

probe_to_gene <- probe_to_gene[
    !is.na(probe_to_gene$symbol) & probe_to_gene$symbol != "",
    ,
    drop = FALSE
]

common_probes <- intersect(
    rownames(expr_raw_gse),
    probe_to_gene$probe_id
)

expr_gse_probe <- expr_raw_gse[
    common_probes,
    ,
    drop = FALSE
]

symbols_gse <- probe_to_gene$symbol[
    match(common_probes, probe_to_gene$probe_id)
]

expr_gse_all <- collapse_probes_by_symbol(
    expr_gse_probe,
    symbols_gse
)

sepsis_ids_gse <- rownames(meta_gse)

expr_gse <- expr_gse_all[
    ,
    sepsis_ids_gse,
    drop = FALSE
]

scores_gse <- score_constructs(expr_gse)

# E-MTAB-5273
emtab_expr_file <- obtain_public_file(
    label = "E-MTAB-5273 processed expression matrix",
    url = paste0(
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-5273/",
        "Burnham_sepsis_discovery_normalised_231.txt"
    ),
    local_dirs = EMTAB_LOCAL_DIRS,
    local_filenames = "Burnham_sepsis_discovery_normalised_231.txt",
    download_file = file.path(
        PROJECT_DIR,
        "Burnham_sepsis_discovery_normalised_231.txt"
    ),
    min_bytes = 10000L
)

emtab_sdrf_file <- obtain_public_file(
    label = "E-MTAB-5273 SDRF metadata",
    url = paste0(
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-5273/",
        "E-MTAB-5273.sdrf.txt"
    ),
    local_dirs = EMTAB_LOCAL_DIRS,
    local_filenames = "E-MTAB-5273.sdrf.txt",
    download_file = file.path(
        PROJECT_DIR,
        "E-MTAB-5273.sdrf.txt"
    ),
    min_bytes = 1000L
)

cat(
    "\nE-MTAB-5273 files selected:\n",
    "  Expression: ", emtab_expr_file, "\n",
    "  SDRF:       ", emtab_sdrf_file, "\n",
    sep = ""
)

raw_expr_emtab <- read.delim(
    emtab_expr_file,
    row.names = 1,
    check.names = FALSE
)

addr_env <- illuminaHumanv4ARRAYADDRESS
probe_ids <- ls(addr_env)

probe_to_addr <- data.frame(
    PROBEID = probe_ids,
    ArrayAddress = vapply(
        probe_ids,
        function(p) as.character(get(p, envir = addr_env)),
        character(1)
    ),
    stringsAsFactors = FALSE
)

addr_matched <- probe_to_addr[
    probe_to_addr$ArrayAddress %in% rownames(raw_expr_emtab),
    ,
    drop = FALSE
]

symbol_map <- AnnotationDbi::select(
    illuminaHumanv4.db,
    keys = unique(addr_matched$PROBEID),
    keytype = "PROBEID",
    columns = "SYMBOL"
)

full_map <- merge(
    addr_matched,
    symbol_map,
    by = "PROBEID"
)

full_map <- full_map[
    !is.na(full_map$SYMBOL) & full_map$SYMBOL != "",
    ,
    drop = FALSE
]

expr_with_symbol <- merge(
    full_map[, c("ArrayAddress", "SYMBOL")],
    data.frame(
        ArrayAddress = rownames(raw_expr_emtab),
        raw_expr_emtab,
        check.names = FALSE
    ),
    by = "ArrayAddress"
)

sample_cols <- setdiff(
    colnames(expr_with_symbol),
    c("ArrayAddress", "SYMBOL")
)

expr_emtab_df <- expr_with_symbol |>
    dplyr::group_by(SYMBOL) |>
    dplyr::summarise(
        dplyr::across(dplyr::all_of(sample_cols), median),
        .groups = "drop"
    ) |>
    as.data.frame()

rownames(expr_emtab_df) <- expr_emtab_df$SYMBOL
expr_emtab_df$SYMBOL <- NULL
expr_emtab <- as.matrix(expr_emtab_df)

sdrf <- read.delim(
    emtab_sdrf_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

normalise_sdrf_name <- function(x) {
    gsub("[^a-z0-9]+", "", tolower(x))
}

find_sdrf_column <- function(df, candidates, label) {
    nm <- colnames(df)

    # 1. Exact match first
    exact <- candidates[candidates %in% nm]
    if (length(exact) > 0L) {
        return(exact[1])
    }

    # 2. Match after removing spaces, punctuation and case differences
    nm_norm <- normalise_sdrf_name(nm)
    cand_norm <- normalise_sdrf_name(candidates)

    idx <- which(nm_norm %in% cand_norm)

    if (length(idx) == 1L) {
        return(nm[idx])
    }

    if (length(idx) > 1L) {
        stop(
            "Multiple SDRF columns matched the requested field '",
            label,
            "': ",
            paste(nm[idx], collapse = ", "),
            "\nPlease inspect colnames(sdrf)."
        )
    }

    stop(
        "Could not identify the SDRF column for '",
        label,
        "'.\nAvailable SDRF columns include:\n  ",
        paste(head(nm, 30), collapse = "\n  "),
        if (length(nm) > 30L) "\n  ..." else ""
    )
}

source_col <- find_sdrf_column(
    sdrf,
    candidates = c(
        "Source Name",
        "Source.Name",
        "SourceName"
    ),
    label = "sample/source name"
)

disease_col <- find_sdrf_column(
    sdrf,
    candidates = c(
        "FactorValue[disease]",
        "Factor Value[disease]",
        "FactorValue.disease.",
        "Factor.Value.disease.",
        "FactorValue disease",
        "disease"
    ),
    label = "disease"
)

cat(
    "\nE-MTAB-5273 SDRF columns selected:\n",
    "  Sample/source column: ", source_col, "\n",
    "  Disease column:       ", disease_col, "\n",
    sep = ""
)

meta_emtab <- data.frame(
    sample_id = as.character(sdrf[[source_col]]),
    disease = as.character(sdrf[[disease_col]]),
    stringsAsFactors = FALSE
)

if (nrow(meta_emtab) == 0L || all(is.na(meta_emtab$sample_id))) {
    stop("E-MTAB-5273 SDRF parsing produced no valid sample identifiers.")
}

meta_emtab$sample_id <- trimws(meta_emtab$sample_id)
meta_emtab$disease <- trimws(meta_emtab$disease)

meta_emtab$patient_id <- stringr::str_extract(
    meta_emtab$sample_id,
    "^[A-Z]+\\d+"
)

time_part <- stringr::str_extract(
    meta_emtab$sample_id,
    "\\.(\\d)$"
)

meta_emtab$timepoint <- ifelse(
    is.na(time_part),
    "1",
    stringr::str_remove(time_part, "\\.")
)

if (anyNA(meta_emtab$patient_id)) {
    bad_ids <- unique(meta_emtab$sample_id[is.na(meta_emtab$patient_id)])
    stop(
        "Could not derive patient_id from one or more E-MTAB-5273 sample names.\n",
        "Examples: ",
        paste(head(bad_ids, 10), collapse = ", ")
    )
}

meta_emtab <- meta_emtab[
    meta_emtab$sample_id %in% colnames(expr_emtab),
    ,
    drop = FALSE
]

meta_emtab <- meta_emtab[
    !is.na(meta_emtab$disease) &
        tolower(meta_emtab$disease) != "normal",
    ,
    drop = FALSE
]

if (nrow(meta_emtab) == 0L) {
    stop(
        "No non-normal E-MTAB-5273 samples remained after matching SDRF ",
        "sample identifiers to the expression-matrix columns."
    )
}

required_meta_cols <- c("sample_id", "timepoint", "patient_id", "disease")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_emtab))

if (length(missing_meta_cols) > 0L) {
    stop(
        "Internal E-MTAB-5273 metadata construction error. Missing column(s): ",
        paste(missing_meta_cols, collapse = ", ")
    )
}

cat(
    sprintf(
        "E-MTAB-5273 metadata prepared: %d non-normal samples, %d patients.\n",
        nrow(meta_emtab),
        length(unique(meta_emtab$patient_id))
    )
)

cat(
    "E-MTAB-5273 nominal timepoints after matching: ",
    paste(
        names(sort(table(meta_emtab$timepoint))),
        as.integer(sort(table(meta_emtab$timepoint))),
        sep = "=",
        collapse = ", "
    ),
    "\n",
    sep = ""
)

dup_emtab <- table(meta_emtab$patient_id, meta_emtab$timepoint)
if (any(dup_emtab > 1L)) {
    stop(
        "E-MTAB-5273 contains duplicated patient x timepoint entries. ",
        "Manual adjudication is required before analysis."
    )
}

scores_emtab <- score_constructs(expr_emtab)

# Analysis datasets
# GSE236713: D1, D2, D5
gse_long <- merge(
    meta_gse[
        meta_gse$timepoint %in% c("D1", "D2", "D5"),
        c("geo_accession", "timepoint", "patient_id"),
        drop = FALSE
    ],
    scores_gse,
    by.x = "geo_accession",
    by.y = "Sample"
)

gse_long$time_num <- as.integer(
    sub("D", "", gse_long$timepoint)
)

# E-MTAB-5273: D1, D3, D5
emtab_long <- merge(
    meta_emtab[
        meta_emtab$timepoint %in% c("1", "3", "5"),
        c("sample_id", "timepoint", "patient_id"),
        drop = FALSE
    ],
    scores_emtab,
    by.x = "sample_id",
    by.y = "Sample"
)

emtab_long$time_num <- as.integer(emtab_long$timepoint)

gse_d1 <- gse_long[
    gse_long$time_num == 1,
    ,
    drop = FALSE
]

emtab_d1 <- emtab_long[
    emtab_long$time_num == 1,
    ,
    drop = FALSE
]

construct_genes <- unique(c(hla_19, tcell_14))

gse_d1_expression <- expr_gse[
    construct_genes,
    gse_d1$geo_accession,
    drop = FALSE
]

emtab_d1_expression <- expr_emtab[
    construct_genes,
    emtab_d1$sample_id,
    drop = FALSE
]

stopifnot(
    ncol(gse_d1_expression) == nrow(gse_d1),
    ncol(emtab_d1_expression) == nrow(emtab_d1)
)

cat("\nSample summary:\n")
cat(
    sprintf(
        "  E-MTAB-5273: D1=%d, D3=%d, D5=%d; unique patients=%d\n",
        sum(emtab_long$time_num == 1),
        sum(emtab_long$time_num == 3),
        sum(emtab_long$time_num == 5),
        length(unique(emtab_long$patient_id))
    )
)
cat(
    sprintf(
        "  GSE236713:   D1=%d, D2=%d, D5=%d; unique patients=%d\n",
        sum(gse_long$time_num == 1),
        sum(gse_long$time_num == 2),
        sum(gse_long$time_num == 5),
        length(unique(gse_long$patient_id))
    )
)

# Raw associations
make_raw_crosssection_row <- function(
    x,
    y,
    cohort,
    role,
    timepoint,
    boot_seed
) {
    ci <- bootstrap_spearman_ci(
        x,
        y,
        B = N_BOOT,
        seed = boot_seed
    )

    tst <- spearman_test(x, y)

    data.frame(
        cohort = cohort,
        role = role,
        timepoint = timepoint,
        n = tst$n,
        rho = unname(ci["rho"]),
        CI_low = unname(ci["lower"]),
        CI_high = unname(ci["upper"]),
        p_value = tst$p,
        stringsAsFactors = FALSE
    )
}

crosssec_at_tp <- function(df, tp, cohort, role, boot_seed) {
    sub <- df[df$time_num == tp, , drop = FALSE]

    make_raw_crosssection_row(
        sub$Score_HLA,
        sub$Score_Tcell,
        cohort = cohort,
        role = role,
        timepoint = paste0("D", tp),
        boot_seed = boot_seed
    )
}

fixed_timepoint_raw <- rbind(
    crosssec_at_tp(
        emtab_long, 1, "E-MTAB-5273", "Primary",
        FIXED_BOOT_SEEDS[["E-MTAB-5273_D1"]]
    ),
    crosssec_at_tp(
        emtab_long, 3, "E-MTAB-5273", "Primary",
        FIXED_BOOT_SEEDS[["E-MTAB-5273_D3"]]
    ),
    crosssec_at_tp(
        emtab_long, 5, "E-MTAB-5273", "Primary",
        FIXED_BOOT_SEEDS[["E-MTAB-5273_D5"]]
    ),
    crosssec_at_tp(
        gse_long, 1, "GSE236713", "Validation",
        FIXED_BOOT_SEEDS[["GSE236713_D1"]]
    ),
    crosssec_at_tp(
        gse_long, 2, "GSE236713", "Validation",
        FIXED_BOOT_SEEDS[["GSE236713_D2"]]
    ),
    crosssec_at_tp(
        gse_long, 5, "GSE236713", "Validation",
        FIXED_BOOT_SEEDS[["GSE236713_D5"]]
    )
)

raw_d1 <- fixed_timepoint_raw[
    fixed_timepoint_raw$timepoint == "D1",
    ,
    drop = FALSE
]

write.csv(
    raw_d1,
    file.path(OUTPUT_DIR, "D1_raw_crosssection.csv"),
    row.names = FALSE
)

write.csv(
    fixed_timepoint_raw,
    file.path(
        OUTPUT_DIR,
        "supportive_fixed_timepoint_raw_correlations.csv"
    ),
    row.names = FALSE
)

cat("\nRaw cross-sectional analyses complete.\n")

# MCPcounter estimates
mcp_clean_gse <- mcp_estimate_local(
    expr_gse,
    MCP_MARKERS_CLEAN
)

mcp_clean_emtab <- mcp_estimate_local(
    expr_emtab,
    MCP_MARKERS_CLEAN
)

mcp_original_gse <- mcp_estimate_local(
    expr_gse,
    MCP_MARKERS_IMMUNE8
)

mcp_original_emtab <- mcp_estimate_local(
    expr_emtab,
    MCP_MARKERS_IMMUNE8
)

stopifnot(
    all(IMMUNE8 %in% rownames(mcp_clean_gse)),
    all(IMMUNE8 %in% rownames(mcp_clean_emtab)),
    all(IMMUNE8 %in% rownames(mcp_original_gse)),
    all(IMMUNE8 %in% rownames(mcp_original_emtab))
)

# Paired changes
paired_changes <- list(
    make_paired_changes(
        emtab_long, "sample_id", 3,
        cohort = "E-MTAB-5273", interval = "D1-D3"
    ),
    make_paired_changes(
        emtab_long, "sample_id", 5,
        cohort = "E-MTAB-5273", interval = "D1-D5"
    ),
    make_paired_changes(
        gse_long, "geo_accession", 2,
        cohort = "GSE236713", interval = "D1-D2"
    ),
    make_paired_changes(
        gse_long, "geo_accession", 5,
        cohort = "GSE236713", interval = "D1-D5"
    )
)

paired_change_data <- do.call(rbind, paired_changes)
rownames(paired_change_data) <- NULL

raw_delta_rows <- vector("list", length(paired_changes))

for (i in seq_along(paired_changes)) {
    pair <- paired_changes[[i]]
    key <- paste(pair$cohort[1], pair$interval[1], sep = "_")

    ci <- bootstrap_spearman_ci(
        pair$dHLA,
        pair$dTcell,
        B = N_BOOT,
        seed = DELTA_BOOT_SEEDS[[key]]
    )

    tst <- spearman_test(pair$dHLA, pair$dTcell)

    raw_delta_rows[[i]] <- data.frame(
        cohort = pair$cohort[1],
        interval = pair$interval[1],
        n = nrow(pair),
        rho = unname(ci["rho"]),
        CI_low = unname(ci["lower"]),
        CI_high = unname(ci["upper"]),
        p_value = tst$p,
        stringsAsFactors = FALSE
    )
}

raw_delta <- do.call(rbind, raw_delta_rows)

write.csv(
    raw_delta,
    file.path(OUTPUT_DIR, "paired_delta_raw.csv"),
    row.names = FALSE
)

# Composition-adjusted analyses
emtab_d1_cov_clean <- get_mcp_covariates(
    mcp_clean_emtab,
    emtab_d1$sample_id
)

gse_d1_cov_clean <- get_mcp_covariates(
    mcp_clean_gse,
    gse_d1$geo_accession
)

adj_primary_rows <- list(
    cbind(
        data.frame(
            cohort = "E-MTAB-5273",
            interval = "D1",
            specification = "Immune8 (clean markers; primary)",
            stringsAsFactors = FALSE
        ),
        run_adjusted_analysis(
            emtab_d1$Score_HLA,
            emtab_d1$Score_Tcell,
            emtab_d1_cov_clean,
            boot_seed = ADJUSTED_BOOT_SEEDS[["E-MTAB-5273_D1"]],
            perm_seed = ADJUSTED_PERM_SEEDS[["E-MTAB-5273_D1"]]
        )
    ),
    cbind(
        data.frame(
            cohort = "GSE236713",
            interval = "D1",
            specification = "Immune8 (clean markers; primary)",
            stringsAsFactors = FALSE
        ),
        run_adjusted_analysis(
            gse_d1$Score_HLA,
            gse_d1$Score_Tcell,
            gse_d1_cov_clean,
            boot_seed = ADJUSTED_BOOT_SEEDS[["GSE236713_D1"]],
            perm_seed = ADJUSTED_PERM_SEEDS[["GSE236713_D1"]]
        )
    )
)

clean_mcp_matrices <- list(
    mcp_clean_emtab,
    mcp_clean_emtab,
    mcp_clean_gse,
    mcp_clean_gse
)

for (i in seq_along(paired_changes)) {
    pair <- paired_changes[[i]]
    key <- paste(pair$cohort[1], pair$interval[1], sep = "_")
    delta_mcp <- make_delta_mcp(pair, clean_mcp_matrices[[i]])

    adj_primary_rows[[length(adj_primary_rows) + 1L]] <- cbind(
        data.frame(
            cohort = pair$cohort[1],
            interval = pair$interval[1],
            specification = "Immune8 (clean markers; primary)",
            stringsAsFactors = FALSE
        ),
        run_adjusted_analysis(
            pair$dHLA,
            pair$dTcell,
            delta_mcp,
            boot_seed = ADJUSTED_BOOT_SEEDS[[key]],
            perm_seed = ADJUSTED_PERM_SEEDS[[key]]
        )
    )
}

adjusted_primary <- do.call(rbind, adj_primary_rows)

desired_order <- c(
    "E-MTAB-5273__D1",
    "E-MTAB-5273__D1-D3",
    "E-MTAB-5273__D1-D5",
    "GSE236713__D1",
    "GSE236713__D1-D2",
    "GSE236713__D1-D5"
)

adjusted_primary$.key <- paste(
    adjusted_primary$cohort,
    adjusted_primary$interval,
    sep = "__"
)

adjusted_primary <- adjusted_primary[
    match(desired_order, adjusted_primary$.key),
    ,
    drop = FALSE
]

adjusted_primary$.key <- NULL

write.csv(
    adjusted_primary,
    file.path(OUTPUT_DIR, "composition_adjusted_primary_clean_Immune8.csv"),
    row.names = FALSE
)

# Marker-set sensitivity
emtab_d1_cov_original <- get_mcp_covariates(
    mcp_original_emtab,
    emtab_d1$sample_id
)

gse_d1_cov_original <- get_mcp_covariates(
    mcp_original_gse,
    gse_d1$geo_accession
)

sens_rows <- list(
    cbind(
        data.frame(
            cohort = "E-MTAB-5273",
            interval = "D1",
            specification = "Immune8 (original markers)",
            stringsAsFactors = FALSE
        ),
        run_sensitivity_analysis(
            emtab_d1$Score_HLA,
            emtab_d1$Score_Tcell,
            emtab_d1_cov_original,
            boot_seed = SENSITIVITY_BOOT_SEEDS[["E-MTAB-5273_D1"]]
        )
    ),
    cbind(
        data.frame(
            cohort = "GSE236713",
            interval = "D1",
            specification = "Immune8 (original markers)",
            stringsAsFactors = FALSE
        ),
        run_sensitivity_analysis(
            gse_d1$Score_HLA,
            gse_d1$Score_Tcell,
            gse_d1_cov_original,
            boot_seed = SENSITIVITY_BOOT_SEEDS[["GSE236713_D1"]]
        )
    )
)

original_mcp_matrices <- list(
    mcp_original_emtab,
    mcp_original_emtab,
    mcp_original_gse,
    mcp_original_gse
)

for (i in seq_along(paired_changes)) {
    pair <- paired_changes[[i]]
    key <- paste(pair$cohort[1], pair$interval[1], sep = "_")
    delta_mcp <- make_delta_mcp(pair, original_mcp_matrices[[i]])

    sens_rows[[length(sens_rows) + 1L]] <- cbind(
        data.frame(
            cohort = pair$cohort[1],
            interval = pair$interval[1],
            specification = "Immune8 (original markers)",
            stringsAsFactors = FALSE
        ),
        run_sensitivity_analysis(
            pair$dHLA,
            pair$dTcell,
            delta_mcp,
            boot_seed = SENSITIVITY_BOOT_SEEDS[[key]]
        )
    )
}

adjusted_original <- do.call(rbind, sens_rows)

adjusted_original$.key <- paste(
    adjusted_original$cohort,
    adjusted_original$interval,
    sep = "__"
)

adjusted_original <- adjusted_original[
    match(desired_order, adjusted_original$.key),
    ,
    drop = FALSE
]

adjusted_original$.key <- NULL

write.csv(
    adjusted_original,
    file.path(OUTPUT_DIR, "composition_adjusted_sensitivity_original_Immune8.csv"),
    row.names = FALSE
)

# Linear mixed-effects models
m_hla_emtab <- lmerTest::lmer(
    Score_HLA ~ time_num + (1 | patient_id),
    data = emtab_long,
    REML = TRUE
)

m_tcell_emtab <- lmerTest::lmer(
    Score_Tcell ~ time_num + (1 | patient_id),
    data = emtab_long,
    REML = TRUE
)

m_hla_gse <- lmerTest::lmer(
    Score_HLA ~ time_num + (1 | patient_id),
    data = gse_long,
    REML = TRUE
)

m_tcell_gse <- lmerTest::lmer(
    Score_Tcell ~ time_num + (1 | patient_id),
    data = gse_long,
    REML = TRUE
)

extract_lmm <- function(model, cohort, role, construct) {
    sm <- summary(model)$coefficients

    data.frame(
        cohort = cohort,
        role = role,
        construct = construct,
        slope_per_day = sm["time_num", "Estimate"],
        SE = sm["time_num", "Std. Error"],
        df = if ("df" %in% colnames(sm)) sm["time_num", "df"] else NA_real_,
        p_value = if ("Pr(>|t|)" %in% colnames(sm)) {
            sm["time_num", "Pr(>|t|)"]
        } else {
            NA_real_
        },
        n_observations = nobs(model),
        n_patients = length(unique(model@frame$patient_id)),
        stringsAsFactors = FALSE
    )
}

lmm_summary <- rbind(
    extract_lmm(
        m_hla_emtab,
        "E-MTAB-5273",
        "Primary",
        "HLA-II antigen-presentation transcriptional programme"
    ),
    extract_lmm(
        m_tcell_emtab,
        "E-MTAB-5273",
        "Primary",
        "T-cell dysfunction-associated transcriptional signature"
    ),
    extract_lmm(
        m_hla_gse,
        "GSE236713",
        "Validation",
        "HLA-II antigen-presentation transcriptional programme"
    ),
    extract_lmm(
        m_tcell_gse,
        "GSE236713",
        "Validation",
        "T-cell dysfunction-associated transcriptional signature"
    )
)

write.csv(
    lmm_summary,
    file.path(OUTPUT_DIR, "supportive_LMM_slopes.csv"),
    row.names = FALSE
)

extract_trajectory <- function(model, days, cohort, role, construct) {
    emm <- emmeans::emmeans(
        model,
        specs = ~ time_num,
        at = list(time_num = days)
    )

    d <- as.data.frame(emm)

    data.frame(
        cohort = cohort,
        role = role,
        construct = construct,
        day = d$time_num,
        timepoint = paste0("D", d$time_num),
        estimated_score = d$emmean,
        SE = d$SE,
        df = d$df,
        CI_low = d$lower.CL,
        CI_high = d$upper.CL,
        stringsAsFactors = FALSE
    )
}

lmm_trajectories <- rbind(
    extract_trajectory(
        m_hla_emtab,
        c(1, 3, 5),
        "E-MTAB-5273",
        "Primary",
        "HLA-II antigen-presentation transcriptional programme"
    ),
    extract_trajectory(
        m_tcell_emtab,
        c(1, 3, 5),
        "E-MTAB-5273",
        "Primary",
        "T-cell dysfunction-associated transcriptional signature"
    ),
    extract_trajectory(
        m_hla_gse,
        c(1, 2, 5),
        "GSE236713",
        "Validation",
        "HLA-II antigen-presentation transcriptional programme"
    ),
    extract_trajectory(
        m_tcell_gse,
        c(1, 2, 5),
        "GSE236713",
        "Validation",
        "T-cell dysfunction-associated transcriptional signature"
    )
)

write.csv(
    lmm_trajectories,
    file.path(OUTPUT_DIR, "supportive_LMM_trajectories.csv"),
    row.names = FALSE
)

# Combined results
ROLE_MAP <- c(
    "E-MTAB-5273" = "Primary",
    "GSE236713" = "Validation"
)

master_results <- rbind(
    # A. Six fixed-timepoint cross-sectional associations
    data.frame(
        section =
            "A. Unadjusted fixed-timepoint cross-sectional associations",
        cohort = fixed_timepoint_raw$cohort,
        role = unname(ROLE_MAP[fixed_timepoint_raw$cohort]),
        interval = fixed_timepoint_raw$timepoint,
        analysis =
            "HLA-II programme vs T-cell dysfunction-associated signature",
        specification = "Unadjusted",
        effect_measure = "Spearman rho",
        n = fixed_timepoint_raw$n,
        n_observations = NA_integer_,
        k = NA_integer_,
        estimate = fixed_timepoint_raw$rho,
        SE = NA_real_,
        CI_low = fixed_timepoint_raw$CI_low,
        CI_high = fixed_timepoint_raw$CI_high,
        p_value = fixed_timepoint_raw$p_value,
        p_method = "Spearman",
        stringsAsFactors = FALSE
    ),

    # B. Four within-patient paired-change associations
    data.frame(
        section =
            "B. Unadjusted within-patient change associations",
        cohort = raw_delta$cohort,
        role = unname(ROLE_MAP[raw_delta$cohort]),
        interval = raw_delta$interval,
        analysis =
            "Delta HLA-II programme vs Delta T-cell dysfunction-associated signature",
        specification = "Unadjusted",
        effect_measure = "Spearman rho",
        n = raw_delta$n,
        n_observations = NA_integer_,
        k = NA_integer_,
        estimate = raw_delta$rho,
        SE = NA_real_,
        CI_low = raw_delta$CI_low,
        CI_high = raw_delta$CI_high,
        p_value = raw_delta$p_value,
        p_method = "Spearman",
        stringsAsFactors = FALSE
    ),

    # C. Six primary clean-marker Immune8-adjusted associations
    data.frame(
        section =
            "C. Primary composition-adjusted associations",
        cohort = adjusted_primary$cohort,
        role = unname(ROLE_MAP[adjusted_primary$cohort]),
        interval = adjusted_primary$interval,
        analysis =
            "HLA-II programme vs T-cell dysfunction-associated signature",
        specification = adjusted_primary$specification,
        effect_measure = "Partial Spearman rho",
        n = adjusted_primary$n,
        n_observations = NA_integer_,
        k = adjusted_primary$k,
        estimate = adjusted_primary$rho,
        SE = NA_real_,
        CI_low = adjusted_primary$CI_low,
        CI_high = adjusted_primary$CI_high,
        p_value = adjusted_primary$p_value,
        p_method = "Freedman-Lane permutation",
        stringsAsFactors = FALSE
    ),

    # D. Six original-marker Immune8 sensitivity rows
    data.frame(
        section =
            "D. Marker-cleaning sensitivity",
        cohort = adjusted_original$cohort,
        role = unname(ROLE_MAP[adjusted_original$cohort]),
        interval = adjusted_original$interval,
        analysis =
            "HLA-II programme vs T-cell dysfunction-associated signature",
        specification = adjusted_original$specification,
        effect_measure = "Partial Spearman rho",
        n = adjusted_original$n,
        n_observations = NA_integer_,
        k = adjusted_original$k,
        estimate = adjusted_original$rho,
        SE = NA_real_,
        CI_low = adjusted_original$CI_low,
        CI_high = adjusted_original$CI_high,
        p_value = NA_real_,
        p_method = "Not reported for targeted sensitivity",
        stringsAsFactors = FALSE
    ),

    # E. Four supportive linear mixed-effects models
    data.frame(
        section =
            "E. Supportive linear mixed-effects models",
        cohort = lmm_summary$cohort,
        role = unname(ROLE_MAP[lmm_summary$cohort]),
        interval = NA_character_,
        analysis = lmm_summary$construct,
        specification = "Score ~ time + (1 | patient)",
        effect_measure = "Beta per day",
        n = lmm_summary$n_patients,
        n_observations = lmm_summary$n_observations,
        k = NA_integer_,
        estimate = lmm_summary$slope_per_day,
        SE = lmm_summary$SE,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_value = lmm_summary$p_value,
        p_method = "Satterthwaite approximation",
        stringsAsFactors = FALSE
    )
)

expected_section_counts <- c(
    "A. Unadjusted fixed-timepoint cross-sectional associations" = 6L,
    "B. Unadjusted within-patient change associations" = 4L,
    "C. Primary composition-adjusted associations" = 6L,
    "D. Marker-cleaning sensitivity" = 6L,
    "E. Supportive linear mixed-effects models" = 4L
)

observed_section_counts <- table(master_results$section)

for (nm in names(expected_section_counts)) {
    observed_n <- if (nm %in% names(observed_section_counts)) {
        as.integer(observed_section_counts[[nm]])
    } else {
        0L
    }

    if (observed_n != expected_section_counts[[nm]]) {
        stop(
            "Master-results section count mismatch for '",
            nm,
            "': expected ",
            expected_section_counts[[nm]],
            ", observed ",
            observed_n,
            "."
        )
    }
}

if (nrow(master_results) != 26L) {
    stop(
        "Expected 26 rows in analysis_master_results.csv; observed ",
        nrow(master_results),
        "."
    )
}

write.csv(
    master_results,
    file.path(
        OUTPUT_DIR,
        "analysis_master_results.csv"
    ),
    row.names = FALSE
)

cat("\nCombined result section counts:\n")
print(observed_section_counts)

cat(
    "\nPASS: analysis_master_results.csv contains the expected 26 rows.\n",
    sep = ""
)

# Session information
analysis_config <- data.frame(
    item = c(
        "Primary cohort",
        "Validation cohort",
        "HLA-II construct",
        "T-cell construct",
        "Primary composition model",
        "Composition sensitivity",
        "Bootstrap replicates",
        "Permutation replicates",
        "GSE236713 preprocessing",
        "Data-access policy"
    ),
    value = c(
        "E-MTAB-5273",
        "GSE236713",
        "19-gene HLA-II antigen-presentation transcriptional programme",
        "14-gene T-cell dysfunction-associated transcriptional signature",
        "MCPcounter Immune8, markers overlapping study constructs removed",
        "MCPcounter Immune8 with original markers",
        as.character(N_BOOT),
        as.character(PARTIAL_PERM_N),
        "Deposited processed/normalised Series Matrix; no additional between-sample normalisation",
        "Local files first; official repository download if absent"
    ),
    stringsAsFactors = FALSE
)

write.csv(
    analysis_config,
    file.path(OUTPUT_DIR, "analysis_configuration.csv"),
    row.names = FALSE
)

writeLines(
    capture.output(sessionInfo()),
    con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

saveRDS(
    list(
        hla_19 = hla_19,
        tcell_14 = tcell_14,
        MCP_MARKERS_IMMUNE8 = MCP_MARKERS_IMMUNE8,
        MCP_MARKERS_CLEAN = MCP_MARKERS_CLEAN,
        IMMUNE8 = IMMUNE8,
        emtab_long = emtab_long,
        gse_long = gse_long,
        emtab_d1 = emtab_d1,
        gse_d1 = gse_d1,
        emtab_d1_expression = emtab_d1_expression,
        gse_d1_expression = gse_d1_expression,
        paired_change_data = paired_change_data,
        sdrf = sdrf,
        meta_gse = meta_gse,
        raw_d1 = raw_d1,
        fixed_timepoint_raw = fixed_timepoint_raw,
        raw_delta = raw_delta,
        adjusted_primary = adjusted_primary,
        adjusted_original = adjusted_original,
        lmm_summary = lmm_summary,
        lmm_trajectories = lmm_trajectories,
        master_results = master_results
    ),
    file = file.path(OUTPUT_DIR, "analysis.rds")
)

cat("\nAnalysis complete.\n")
cat(
    "Outputs written to:\n  ",
    normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE),
    "\n",
    sep = ""
)
