analysis_file <- file.path("outputs_publication", "analysis.rds")

if (!file.exists(analysis_file)) {
    stop("Required analysis output was not found. Run analysis.R or run.R first.")
}

obj <- readRDS(analysis_file)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

hla_genes <- obj$hla_19
tcell_genes <- obj$tcell_14
emtab_expr <- obj$emtab_d1_expression
gse_expr <- obj$gse_d1_expression

stopifnot(length(hla_genes) == 19L, length(tcell_genes) == 14L)

make_rows <- function(genes, construct) {
    data.frame(
        Construct = construct,
        Gene = genes,
        `E-MTAB-5273` = ifelse(
            genes %in% rownames(emtab_expr),
            "Available",
            "Not available"
        ),
        GSE236713 = ifelse(
            genes %in% rownames(gse_expr),
            "Available",
            "Not available"
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

table_s2 <- rbind(
    make_rows(
        hla_genes,
        "HLA-II antigen-presentation transcriptional programme (19 genes)"
    ),
    make_rows(
        tcell_genes,
        "T-cell dysfunction-associated transcriptional signature (14 genes)"
    )
)

availability_summary <- data.frame(
    Cohort = c("E-MTAB-5273", "GSE236713"),
    `HLA-II genes available` = c(
        sum(hla_genes %in% rownames(emtab_expr)),
        sum(hla_genes %in% rownames(gse_expr))
    ),
    `HLA-II total` = c(length(hla_genes), length(hla_genes)),
    `T-cell signature genes available` = c(
        sum(tcell_genes %in% rownames(emtab_expr)),
        sum(tcell_genes %in% rownames(gse_expr))
    ),
    `T-cell signature total` = c(length(tcell_genes), length(tcell_genes)),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (any(table_s2$`E-MTAB-5273` != "Available") ||
    any(table_s2$GSE236713 != "Available")) {
    stop("One or more study genes are unavailable in a D1 expression matrix.")
}

write.csv(
    table_s2,
    file.path("tables", "TableS2.csv"),
    row.names = FALSE,
    na = ""
)

write.csv(
    availability_summary,
    file.path("tables", "TableS2_availability_summary.csv"),
    row.names = FALSE,
    na = ""
)

print(table_s2, row.names = FALSE)
