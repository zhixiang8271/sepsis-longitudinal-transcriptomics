analysis_file <- file.path("outputs_publication", "analysis.rds")

if (!file.exists(analysis_file)) {
    stop("Required analysis output was not found. Run analysis.R or run.R first.")
}

obj <- readRDS(analysis_file)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

normalise_name <- function(x) {
    gsub("[^a-z0-9]+", "", tolower(x))
}

find_column <- function(df, candidates, label) {
    nm <- colnames(df)
    nm_norm <- normalise_name(nm)
    cand_norm <- normalise_name(candidates)
    idx <- which(nm_norm %in% cand_norm)

    if (length(idx) != 1L) {
        stop(
            "Could not identify a unique column for ", label,
            ". Matched: ",
            if (length(idx) == 0L) "none" else paste(nm[idx], collapse = ", ")
        )
    }

    nm[idx]
}

fmt_n_pct <- function(n, denominator) {
    if (is.na(n) || is.na(denominator) || denominator <= 0) {
        return("Not available")
    }

    sprintf("%d (%.1f%%)", n, 100 * n / denominator)
}

fmt_med_iqr <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return("Not available")
    }

    q <- quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
    sprintf("%.1f [%.1f-%.1f]", q[2], q[1], q[3])
}

sdrf <- obj$sdrf
emtab_d1 <- obj$emtab_d1
gse_d1 <- obj$gse_d1
meta_gse <- obj$meta_gse

source_col <- find_column(
    sdrf,
    c("Source Name", "Source.Name", "SourceName"),
    "E-MTAB sample identifier"
)

sex_col <- find_column(
    sdrf,
    c("Characteristics[sex]", "Characteristics.sex.", "sex"),
    "E-MTAB sex"
)

age_col <- find_column(
    sdrf,
    c("Characteristics[age]", "Characteristics.age.", "age"),
    "E-MTAB age"
)

clinical_col <- find_column(
    sdrf,
    c(
        "Characteristics[clinical information]",
        "Characteristics.clinical.information.",
        "clinical information"
    ),
    "E-MTAB clinical information"
)

match_emtab <- match(emtab_d1$sample_id, as.character(sdrf[[source_col]]))
if (anyNA(match_emtab)) {
    stop("One or more E-MTAB D1 samples could not be matched to the SDRF.")
}

sdrf_d1 <- sdrf[match_emtab, , drop = FALSE]
n_emtab <- nrow(emtab_d1)

age_emtab <- suppressWarnings(
    as.numeric(trimws(as.character(sdrf_d1[[age_col]])))
)

sex_raw <- tolower(trimws(as.character(sdrf_d1[[sex_col]])))
sex_emtab <- ifelse(
    sex_raw %in% c("female", "f"),
    "Female",
    ifelse(sex_raw %in% c("male", "m"), "Male", NA_character_)
)

female_emtab <- sum(sex_emtab == "Female", na.rm = TRUE)
male_emtab <- sum(sex_emtab == "Male", na.rm = TRUE)
sex_known_emtab <- sum(!is.na(sex_emtab))

source_emtab <- ifelse(
    grepl("^CAP", emtab_d1$sample_id),
    "CAP",
    ifelse(grepl("^FP", emtab_d1$sample_id), "FP", NA_character_)
)

if (anyNA(source_emtab)) {
    stop("One or more E-MTAB D1 samples could not be classified as CAP or FP.")
}

cap_emtab <- sum(source_emtab == "CAP")
fp_emtab <- sum(source_emtab == "FP")

clinical_raw <- tolower(trimws(as.character(sdrf_d1[[clinical_col]])))
alive_idx <- !is.na(clinical_raw) & grepl("alive.*28|28.*alive", clinical_raw)
dead_idx <- !is.na(clinical_raw) & grepl("(dead|died).*28|28.*(dead|died)", clinical_raw)

if (any(alive_idx & dead_idx, na.rm = TRUE)) {
    stop("An E-MTAB clinical-information entry matched both survival patterns.")
}

outcome_emtab <- rep(NA_character_, n_emtab)
outcome_emtab[alive_idx] <- "Survived"
outcome_emtab[dead_idx] <- "Died"

emtab_outcome_known <- sum(!is.na(outcome_emtab))
emtab_survived <- sum(outcome_emtab == "Survived", na.rm = TRUE)
emtab_died <- sum(outcome_emtab == "Died", na.rm = TRUE)

match_gse <- match(gse_d1$geo_accession, meta_gse$geo_accession)
if (anyNA(match_gse)) {
    stop("One or more GSE236713 D1 samples could not be matched to metadata.")
}

gse_meta_d1 <- meta_gse[match_gse, , drop = FALSE]
n_gse <- nrow(gse_d1)

source_gse <- as.character(gse_meta_d1$infection_source)
source_gse[!source_gse %in% c("Pulmonary", "Abdominal")] <- NA_character_

pulmonary_gse <- sum(source_gse == "Pulmonary", na.rm = TRUE)
abdominal_gse <- sum(source_gse == "Abdominal", na.rm = TRUE)
source_known_gse <- sum(!is.na(source_gse))

outcome_gse <- as.character(gse_meta_d1$outcome_28d)
outcome_gse[!outcome_gse %in% c("Survived", "Died")] <- NA_character_

gse_outcome_known <- sum(!is.na(outcome_gse))
gse_survived <- sum(outcome_gse == "Survived", na.rm = TRUE)
gse_died <- sum(outcome_gse == "Died", na.rm = TRUE)

stopifnot(
    cap_emtab + fp_emtab == n_emtab,
    male_emtab + female_emtab == sex_known_emtab,
    emtab_survived + emtab_died == emtab_outcome_known,
    pulmonary_gse + abdominal_gse == source_known_gse,
    gse_survived + gse_died == gse_outcome_known
)

table_s1 <- data.frame(
    Characteristic = c(
        "Day-1 analytic samples, n",
        "Age, years, median [IQR]",
        "Female sex, n (%)",
        "28-day mortality, n (%)",
        "Infection source, n (%)"
    ),
    `E-MTAB-5273 (primary)` = c(
        as.character(n_emtab),
        fmt_med_iqr(age_emtab),
        fmt_n_pct(female_emtab, sex_known_emtab),
        fmt_n_pct(emtab_died, emtab_outcome_known),
        sprintf(
            "CAP %d (%.1f%%); FP %d (%.1f%%)",
            cap_emtab, 100 * cap_emtab / n_emtab,
            fp_emtab, 100 * fp_emtab / n_emtab
        )
    ),
    `GSE236713 (validation)` = c(
        as.character(n_gse),
        "Not available",
        "Not available",
        fmt_n_pct(gse_died, gse_outcome_known),
        if (source_known_gse > 0L) {
            sprintf(
                "Pulmonary %d (%.1f%%); Abdominal %d (%.1f%%)",
                pulmonary_gse, 100 * pulmonary_gse / source_known_gse,
                abdominal_gse, 100 * abdominal_gse / source_known_gse
            )
        } else {
            "Not available"
        }
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

table_s1_numeric <- data.frame(
    Cohort = c("E-MTAB-5273", "GSE236713"),
    Role = c("Primary", "Validation"),
    D1_n = c(n_emtab, n_gse),
    Age_available_n = c(sum(is.finite(age_emtab)), NA_integer_),
    Age_median = c(median(age_emtab, na.rm = TRUE), NA_real_),
    Age_Q1 = c(unname(quantile(age_emtab, 0.25, na.rm = TRUE)), NA_real_),
    Age_Q3 = c(unname(quantile(age_emtab, 0.75, na.rm = TRUE)), NA_real_),
    Sex_available_n = c(sex_known_emtab, NA_integer_),
    Male_n = c(male_emtab, NA_integer_),
    Female_n = c(female_emtab, NA_integer_),
    Outcome_28d_available_n = c(emtab_outcome_known, gse_outcome_known),
    Survived_28d_n = c(emtab_survived, gse_survived),
    Died_28d_n = c(emtab_died, gse_died),
    CAP_n = c(cap_emtab, NA_integer_),
    FP_n = c(fp_emtab, NA_integer_),
    Pulmonary_n = c(NA_integer_, pulmonary_gse),
    Abdominal_n = c(NA_integer_, abdominal_gse),
    Infection_source_available_n = c(n_emtab, source_known_gse),
    stringsAsFactors = FALSE
)

write.csv(
    table_s1,
    file.path("tables", "TableS1.csv"),
    row.names = FALSE,
    na = ""
)

write.csv(
    table_s1_numeric,
    file.path("tables", "TableS1_numeric.csv"),
    row.names = FALSE,
    na = ""
)

print(table_s1, row.names = FALSE)
