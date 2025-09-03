#!/usr/bin/env Rscript

# Synergy analysis between two genes using public breast cancer datasets
# Automatically downloads data from TCGA-BRCA (UCSC Xena), METABRIC (cBioPortalData)
# and GSE96058 (GEO). Performs simple statistical analyses: Spearman correlation,
# Fisher's exact test for double-high enrichment, KM curves, and Cox regression
# with interaction term.

## ------------------------- Package handling -------------------------------
options(timeout = 600)

required_pkgs_cran <- c(
  "optparse", "UCSCXenaTools", "GEOquery", "SummarizedExperiment",
  "MultiAssayExperiment", "dplyr", "data.table", "survival",
  "survminer", "ggplot2", "readr", "stringr"
)
required_pkgs_bioc <- c("cBioPortalData", "org.Hs.eg.db")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    message("Installing package: ", pkg)
    if (bioc) {
      if (!suppressWarnings(require("BiocManager", quietly = TRUE))) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

for (p in required_pkgs_cran) install_if_missing(p)
for (p in required_pkgs_bioc) install_if_missing(p, bioc = TRUE)

suppressPackageStartupMessages({
  library(optparse)
  library(UCSCXenaTools)
  library(GEOquery)
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  library(cBioPortalData)
  library(dplyr)
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(org.Hs.eg.db)
})

## ------------------------- Utilities --------------------------------------

safe_readr <- function(path, ...) {
  tryCatch(readr::write_lines("", path, append = TRUE), error = function(e) NULL)
}

write_session_info <- function(outdir) {
  dir.create(file.path(outdir, "log"), showWarnings = FALSE, recursive = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "log", "sessionInfo.txt"))
}

save_data_source <- function(outdir, text) {
  dir.create(file.path(outdir, "log"), showWarnings = FALSE, recursive = TRUE)
  writeLines(text, file.path(outdir, "log", "data_sources.txt"))
}

zscore <- function(x) as.numeric(scale(x))

prep_survival <- function(df) {
  # unify time/status fields, assume time in days
  time_cols <- c("OS_TIME", "OS.time", "OS_MONTHS", "OS_MONTH", "OS_MONTHS.x",
                 "overall.survival.days", "PFI.time")
  status_cols <- c("OS_STATUS", "OS.event", "OS_EVENT", "overall.survival.event", "PFI")
  time_col <- intersect(time_cols, names(df))
  status_col <- intersect(status_cols, names(df))
  if (length(time_col) == 0 || length(status_col) == 0) {
    stop("No survival fields found")
  }
  time <- as.numeric(df[[time_col[1]]])
  status_raw <- as.character(df[[status_col[1]]])
  status <- ifelse(grepl("DECEASED|dead|Death|1", status_raw, ignore.case = TRUE), 1, 0)
  if (grepl("MONTH", time_col[1], ignore.case = TRUE)) {
    time <- time * 30.44
  }
  data.frame(time = time, status = status)
}

# dichotomize by median
binary_median <- function(x) {
  ifelse(x > median(x, na.rm = TRUE), 1, 0)
}

## ------------------------- Analysis functions ----------------------------

do_spearman <- function(df, g1, g2) {
  res <- cor.test(df[[g1]], df[[g2]], method = "spearman")
  data.frame(rho = res$estimate, p = res$p.value)
}

do_fisher <- function(df, g1, g2) {
  g1_hi <- binary_median(df[[g1]])
  g2_hi <- binary_median(df[[g2]])
  tab <- table(g1_hi, g2_hi)
  fisher <- fisher.test(tab)
  exp_hh <- (sum(g1_hi) / length(g1_hi)) * (sum(g2_hi) / length(g2_hi))
  obs_hh <- tab[2,2] / sum(tab)
  data.frame(
    g1_hi0_g2_hi0 = tab[1,1],
    g1_hi1_g2_hi0 = tab[2,1],
    g1_hi0_g2_hi1 = tab[1,2],
    g1_hi1_g2_hi1 = tab[2,2],
    expected_hh = exp_hh,
    observed_hh = obs_hh,
    fisher_p = fisher$p.value
  )
}

do_km_plot <- function(df, g1, g2, outpath) {
  g1_hi <- binary_median(df[[g1]])
  g2_hi <- binary_median(df[[g2]])
  grp <- interaction(g1_hi, g2_hi)
  grp <- factor(grp, levels = c("0.0", "1.0", "0.1", "1.1"), labels = c("LL","HL","LH","HH"))
  surv_obj <- Surv(df$time, df$status)
  fit <- survfit(surv_obj ~ grp)
  g <- ggsurvplot(fit, data = data.frame(df, grp = grp), risk.table = TRUE,
                  legend.title = "Group", legend.labs = levels(grp))
  ggsave(outpath, g$plot, width = 6, height = 5, dpi = 300)
  invisible(fit)
}

do_cox_interaction <- function(df, g1, g2, covars = NULL) {
  form <- as.formula(paste0("Surv(time, status) ~ scale(", g1, ")*scale(", g2, ")"))
  if (!is.null(covars)) {
    cv <- paste(covars, collapse = "+")
    form <- as.formula(paste0("Surv(time, status) ~ scale(", g1, ")*scale(", g2, ") + ", cv))
  }
  fit <- coxph(form, data = df)
  int_term <- paste0("scale(", g1, "):scale(", g2, ")")
  co <- summary(fit)$coefficients[int_term,]
  ci <- confint(fit)[int_term,]
  data.frame(
    term = int_term,
    beta = co["coef"],
    HR = exp(co["coef"]),
    lower95 = exp(ci[1]),
    upper95 = exp(ci[2]),
    p = co["Pr(>|z|)"]
  )
}

## ------------------------- Data loaders ----------------------------------

load_tcga <- function(genes) {
  message("Loading TCGA-BRCA from UCSC Xena")
  host <- "https://gdc.xenahubs.net"
  ds <- UCSCXenaTools::XenaDatasets(host)
  expr_ds <- ds[stringr::str_detect(ds$XenaDatasets, "TCGA.BRCA") &
                  stringr::str_detect(ds$XenaDatasets, "gene"), ]
  expr_ds <- expr_ds[stringr::str_detect(expr_ds$XenaDatasets, "FPKM|TPM|RSEM"), ]
  expr_name <- expr_ds$XenaDatasets[1]
  message("Expression dataset:", expr_name)
  exp <- UCSCXenaTools::XenaGenerate(subset = expr_name, host = host) %>%
    UCSCXenaTools::XenaQuery() %>%
    UCSCXenaTools::XenaDownload(destdir = tempdir(), method = "auto") %>%
    UCSCXenaTools::XenaPrepare()
  exp <- as.data.frame(exp)
  exp <- exp[genes, , drop = FALSE]
  exp <- t(exp)
  exp <- as.data.frame(exp)
  exp$sample <- rownames(exp)
  clinical_ds <- ds[stringr::str_detect(ds$XenaDatasets, "TCGA.BRCA") &
                      stringr::str_detect(ds$XenaDatasets, "clinical|phenotype"), ]
  clin_name <- clinical_ds$XenaDatasets[1]
  message("Clinical dataset:", clin_name)
  clin <- UCSCXenaTools::XenaGenerate(subset = clin_name, host = host) %>%
    UCSCXenaTools::XenaQuery() %>%
    UCSCXenaTools::XenaDownload(destdir = tempdir(), method = "auto") %>%
    UCSCXenaTools::XenaPrepare()
  clin <- as.data.frame(clin)
  clin$sample <- rownames(clin)
  df <- merge(exp, clin, by = "sample")
  # sample type filter: TCGA barcodes positions 14-15 "01" for tumor
  df <- df[substr(df$sample, 14, 15) == "01", ]
  surv <- prep_survival(df)
  df$time <- surv$time
  df$status <- surv$status
  df$sample <- NULL
  df <- df[complete.cases(df[, c(genes, "time", "status")]), ]
  df[genes] <- lapply(df[genes], zscore)
  list(data = df, pam50 = df$PAM50.mRNA)
}

load_metabric <- function(genes) {
  message("Loading METABRIC via cBioPortalData")
  mb <- cBioPortalData::cBioDataPack("brca_metabric", use_cache = TRUE, ask = FALSE)
  exp <- assays(mb$mrna)$mrna
  exp <- exp[genes, ]
  exp <- t(exp)
  exp <- as.data.frame(exp)
  exp$sample <- rownames(exp)
  clin <- as.data.frame(colData(mb))
  clin$sample <- rownames(clin)
  df <- merge(exp, clin, by = "sample")
  surv <- prep_survival(df)
  df$time <- surv$time
  df$status <- surv$status
  df$sample <- NULL
  df <- df[complete.cases(df[, c(genes, "time", "status")]), ]
  df[genes] <- lapply(df[genes], zscore)
  list(data = df, pam50 = df$PAM50)
}

load_gse96058 <- function(genes) {
  message("Loading GSE96058 via GEOquery")
  gset <- GEOquery::getGEO("GSE96058", GSEMatrix = TRUE, getGPL = FALSE)
  gset <- gset[[1]]
  exp <- exprs(gset)
  rn <- rownames(exp)
  if (all(grepl("^ENSG", rn))) {
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = rn, columns = "SYMBOL", keytype = "ENSEMBL")
    map <- map[!duplicated(map$ENSEMBL), ]
    rn2 <- map$SYMBOL[match(rn, map$ENSEMBL)]
    keep <- !is.na(rn2)
    exp <- exp[keep, ]
    rownames(exp) <- rn2[keep]
  }
  if (!all(genes %in% rownames(exp))) {
    stop("Some genes not found in GSE96058")
  }
  exp <- exp[genes, ]
  exp <- t(exp)
  exp <- as.data.frame(exp)
  exp$sample <- rownames(exp)
  clin <- pData(gset)
  clin$sample <- rownames(clin)
  df <- merge(exp, clin, by = "sample")
  surv <- prep_survival(df)
  df$time <- surv$time
  df$status <- surv$status
  df$sample <- NULL
  df <- df[complete.cases(df[, c(genes, "time", "status")]), ]
  df[genes] <- lapply(df[genes], zscore)
  list(data = df, pam50 = df$PAM50)
}

## ------------------------- Main workflow ---------------------------------

option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = "USP25,C1QL4",
              help = "Two genes separated by comma"),
  make_option(c("-d", "--datasets"), type = "character", default = "tcga,metabric,gse96058",
              help = "Datasets to run"),
  make_option(c("-o", "--outdir"), type = "character", default = "./synergy_out"),
  make_option(c("-p", "--pam50-strata"), type = "character", default = "false"),
  make_option(c("-m", "--min-samples"), type = "integer", default = 100),
  make_option(c("-s", "--seed"), type = "integer", default = 2025)
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

genes <- strsplit(opt$genes, ",")[[1]]
if (length(genes) != 2) stop("Exactly two genes required")

runs <- strsplit(opt$datasets, ",")[[1]]
run_pam <- tolower(opt$`pam50-strata`) == "true"

outdir <- opt$outdir
summary_list <- list()
summary_overall <- data.frame()

for (ds in runs) {
  ds <- tolower(ds)
  message("Processing dataset: ", ds)
  res <- tryCatch({
    if (ds == "tcga") load_tcga(genes)
    else if (ds == "metabric") load_metabric(genes)
    else if (ds == "gse96058") load_gse96058(genes)
    else stop("Unknown dataset: ", ds)
  }, error = function(e) {
    message("Failed to load ", ds, ": ", e$message)
    return(NULL)
  })
  if (is.null(res)) next
  df <- res$data
  pam <- res$pam50
  if (nrow(df) < opt$`min-samples`) {
    message("Not enough samples in ", ds)
    next
  }
  ds_dir <- file.path(outdir, ds)
  dir.create(file.path(ds_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ds_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  # Spearman
  spearman <- do_spearman(df, genes[1], genes[2])
  readr::write_csv(spearman, file.path(ds_dir, "tables", "correlation.csv"))

  # Fisher
  fisher <- do_fisher(df, genes[1], genes[2])
  readr::write_csv(fisher, file.path(ds_dir, "tables", "fisher_2x2.csv"))

  # KM plot
  km_fit <- try(do_km_plot(df, genes[1], genes[2], file.path(ds_dir, "plots", "km_4groups.png")))

  # Cox interaction
  covars <- intersect(c("age", "AGE", "AGE_AT_DIAGNOSIS", "stage", "subtype"), names(df))
  cox <- do_cox_interaction(df, genes[1], genes[2], covars)
  cox$covariates <- paste(covars, collapse = ",")
  cox$samples <- nrow(df)
  readr::write_csv(cox, file.path(ds_dir, "tables", "cox_interaction.csv"))

  summary_overall <- rbind(summary_overall, data.frame(dataset = ds,
                                                       spearman_rho = spearman$rho,
                                                       spearman_p = spearman$p,
                                                       fisher_p = fisher$fisher_p,
                                                       cox_hr = cox$HR,
                                                       cox_p = cox$p))
  message(sprintf("%s: Spearman rho=%.3f (p=%.3g); Fisher p=%.3g; Cox HR=%.3f (p=%.3g)",
                  ds, spearman$rho, spearman$p, fisher$fisher_p, cox$HR, cox$p))

  # PAM50 stratification
  if (run_pam && !all(is.na(pam))) {
    pam_col <- pam
    df$pam50 <- pam_col
    subtypes <- na.omit(unique(pam_col))
    for (st in subtypes) {
      message("  PAM50 subtype: ", st)
      df_sub <- df[df$pam50 == st, ]
      if (nrow(df_sub) < opt$`min-samples`) {
        message("   Not enough samples")
        next
      }
      st_dir <- file.path(ds_dir, paste0("pam50_", st))
      dir.create(file.path(st_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(st_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
      spearman <- do_spearman(df_sub, genes[1], genes[2])
      readr::write_csv(spearman, file.path(st_dir, "tables", paste0("correlation_", st, ".csv")))
      fisher <- do_fisher(df_sub, genes[1], genes[2])
      readr::write_csv(fisher, file.path(st_dir, "tables", paste0("fisher_", st, ".csv")))
      try(do_km_plot(df_sub, genes[1], genes[2], file.path(st_dir, "plots", paste0("km_", st, ".png"))))
      cox <- do_cox_interaction(df_sub, genes[1], genes[2], covars)
      cox$covariates <- paste(covars, collapse = ",")
      cox$samples <- nrow(df_sub)
      readr::write_csv(cox, file.path(st_dir, "tables", paste0("cox_", st, ".csv")))
      message(sprintf("   %s: Spearman rho=%.3f (p=%.3g); Fisher p=%.3g; Cox HR=%.3f (p=%.3g)",
                      st, spearman$rho, spearman$p, fisher$fisher_p, cox$HR, cox$p))
    }
  }
}

if (nrow(summary_overall) > 0) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(summary_overall, file.path(outdir, "summary_overall.csv"))
  print(summary_overall)
}

write_session_info(outdir)

message("Analysis finished")

