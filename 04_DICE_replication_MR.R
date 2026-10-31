# =============================================================================
# DICE replication MR analysis
# =============================================================================
# Purpose:
#   This script replicates the MR associations supported by the colocalization
#   analysis using immune-cell eQTL data from the DICE database.
#
# Workflow:
#   1. Read positive cell-gene-trait pairs from the colocalization results.
#   2. Map the sc-eQTL cell labels to the corresponding DICE cell labels.
#   3. Read the matching DICE .vcf.gz file and select the target gene.
#   4. Select genome-wide significant DICE eQTL instruments.
#   5. Infer standard errors from DICE beta estimates and P values.
#   6. Harmonize DICE eQTL instruments with the standardized GWAS outcome.
#   7. Run Wald-ratio and inverse-variance-weighted MR where applicable.
#   8. Save pair-level and combined replication results.
#
# DICE input columns:
#   #CHROM, POS, ID, REF, ALT, Gene, GeneSymbol, Pvalue, Beta, Cell
#
# Standardized GWAS outcome columns:
#   SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval, N
# =============================================================================

source("config.R")
library(TwoSampleMR)
library(data.table)

## ---- Parameters -------------------------------------------------------------
DICE <- c("CD4_NAIVE", "CD8_NAIVE", "MONOCYTES", "B_CELL_NAIVE", "NK", "M2")
DICE_MAP <- c(
  cd4nc = "CD4_NAIVE",
  cd8nc = "CD8_NAIVE",
  monoc = "MONOCYTES",
  bin   = "B_CELL_NAIVE",
  nk    = "NK",
  mononc = "M2"
)

P_THRESHOLD <- 5e-8
PPH4_THRESHOLD <- 0.90

COLOC_RESULTS_FILE <- file.path(RESULTS_DIR, "coloc", "coloc_summary_all.tsv")
DICE_DIR <- DATA_DIR
OUTCOME_DIR <- DATA_DIR
RESULTS_DICE_DIR <- file.path(RESULTS_DIR, "dice_replication")
dir.create(RESULTS_DICE_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- Helpers ----------------------------------------------------------------
read_outcome <- function(trait) {
  path <- file.path(OUTCOME_DIR, paste0(trait, ".txt"))
  if (!file.exists(path)) stop("Outcome file not found: ", path)
  required <- c("SNP", "chr", "pos", "effect_allele", "other_allele",
                "eaf", "beta", "se", "pval", "N")
  dat <- fread(path, sep = "\t")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Outcome file is missing columns: ", paste(missing, collapse = ", "))
  }
  dat
}

read_dice <- function(dice_cell) {
  path <- file.path(DICE_DIR, paste0(dice_cell, ".vcf.gz"))
  if (!file.exists(path)) stop("DICE file not found: ", path)

  # skip locates the VCF header line beginning with #CHROM.
  dat <- fread(path, sep = "\t", skip = "#CHROM", header = TRUE,
               check.names = FALSE)
  if ("#CHROM" %in% names(dat)) setnames(dat, "#CHROM", "CHROM")

  required <- c("CHROM", "POS", "ID", "REF", "ALT", "Gene",
                "GeneSymbol", "Pvalue", "Beta", "Cell")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("DICE file is missing columns: ", paste(missing, collapse = ", "))
  }
  dat
}

prepare_dice_exposure <- function(dice_cell, gene_id) {
  dat <- read_dice(dice_cell)
  dat <- dat[Cell == dice_cell & (Gene == gene_id | GeneSymbol == gene_id)]
  if (nrow(dat) == 0) {
    stop("Gene not found in DICE data: ", gene_id, " (", dice_cell, ")")
  }

  dat <- dat[, .(
    SNP = as.character(ID),
    chr = as.character(CHROM),
    pos = as.numeric(POS),
    effect_allele = toupper(as.character(ALT)),
    other_allele = toupper(as.character(REF)),
    beta = as.numeric(Beta),
    pval = as.numeric(Pvalue)
  )]
  dat <- dat[!duplicated(SNP) & is.finite(beta) & is.finite(pval) &
               pval > 0 & pval <= P_THRESHOLD]
  if (nrow(dat) == 0) return(NULL)

  # DICE provides beta and Pvalue but no SE; derive SE from the two-sided P value.
  p_for_z <- pmax(dat$pval, .Machine$double.xmin)
  z_abs <- qnorm(p_for_z / 2, lower.tail = FALSE)
  dat[, se := abs(beta) / z_abs]
  dat <- dat[is.finite(se) & se > 0]
  if (nrow(dat) == 0) return(NULL)

  format_data(
    dat,
    type = "exposure",
    snp_col = "SNP",
    chr_col = "chr",
    pos_col = "pos",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    phenotype_col = "SNP"
  )
}

run_one_replication <- function(source_cell, gene_id, outcome_trait) {
  dice_cell <- unname(DICE_MAP[[source_cell]])
  if (is.null(dice_cell) || is.na(dice_cell)) {
    message("[skip] No DICE mapping for cell: ", source_cell)
    return(NULL)
  }
  message("\n------ DICE replication: ", source_cell, " -> ", dice_cell,
          ": ", gene_id, " -> ", outcome_trait, " ------")

  exposure_dat <- tryCatch(
    prepare_dice_exposure(dice_cell, gene_id),
    error = function(e) {
      message("  [error] DICE exposure: ", e$message)
      NULL
    }
  )
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
    message("  [skip] no genome-wide significant DICE instrument")
    return(NULL)
  }

  outcome_dat <- tryCatch(
    read_outcome(outcome_trait),
    error = function(e) {
      message("  [error] outcome: ", e$message)
      NULL
    }
  )
  if (is.null(outcome_dat)) return(NULL)

  outcome_dat <- outcome_dat[SNP %in% exposure_dat$SNP]
  if (nrow(outcome_dat) == 0) {
    message("  [skip] no DICE instrument found in outcome")
    return(NULL)
  }
  outcome_dat <- format_data(
    outcome_dat,
    type = "outcome",
    snp_col = "SNP",
    chr_col = "chr",
    pos_col = "pos",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    samplesize_col = "N",
    phenotype_col = "outcome"
  )

  exposure_dat$exposure <- exposure_dat$id.exposure <- paste(source_cell, gene_id, sep = ":")
  outcome_dat$outcome <- outcome_dat$id.outcome <- outcome_trait
  dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
  dat <- dat[mr_keep %in% TRUE]
  if (nrow(dat) == 0) {
    message("  [skip] no SNPs after harmonization")
    return(NULL)
  }

  mr_res <- tryCatch(mr(dat), error = function(e) {
    message("  [error] MR: ", e$message)
    NULL
  })
  if (is.null(mr_res) || nrow(mr_res) == 0) return(NULL)

  mr_res[, `:=`(source_cell = source_cell, dice_cell = dice_cell,
                gene = gene_id, outcome_trait = outcome_trait,
                nsnp_harmonized = nrow(dat))]
  pair_dir <- file.path(RESULTS_DICE_DIR, source_cell, gene_id)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(mr_res, file.path(pair_dir, paste0(outcome_trait, ".MR.txt")), sep = "\t")

  mr_res
}

## ---- Read positive coloc pairs ---------------------------------------------
positive <- fread(COLOC_RESULTS_FILE, sep = "\t")
required_coloc <- c("cell", "gene", "trait", "PP.H4")
missing_coloc <- setdiff(required_coloc, names(positive))
if (length(missing_coloc) > 0) {
  stop("Colocalization results are missing columns: ",
       paste(missing_coloc, collapse = ", "))
}
positive <- unique(positive[PP.H4 >= PPH4_THRESHOLD,
                            .(cell, gene, trait)])
positive <- positive[cell %in% names(DICE_MAP) & trait %in% all_traits]

## ---- Run replication for all positive pairs --------------------------------
all_results <- list()
log_rows <- list()

for (i in seq_len(nrow(positive))) {
  pair <- positive[i]
  result <- run_one_replication(pair$cell, pair$gene, pair$trait)
  if (is.null(result)) next
  all_results[[length(all_results) + 1L]] <- result
  log_rows[[length(log_rows) + 1L]] <- data.table(
    source_cell = pair$cell,
    dice_cell = DICE_MAP[[pair$cell]],
    gene = pair$gene,
    outcome_trait = pair$trait,
    status = "success"
  )
}

if (length(all_results) > 0) {
  fwrite(rbindlist(all_results, fill = TRUE),
         file.path(RESULTS_DICE_DIR, "DICE_replication_all_results.txt"),
         sep = "\t")

  # Five physical frailty components are treated as one testing family.
  physical_frailty_traits <- c(
    "Weight_loss",
    "Exhaustion",
    "Low_Physical_activity",
    "Slow_Walking_speed",
    "Low_Grip_strength"
  )

  dice_all[, fdr_family := fifelse(
    outcome_trait %in% physical_frailty_traits,
    "physical_frailty",
    outcome_trait
  )]

  # Apply BH correction separately for each MR method and testing family.
  dice_all[, FDR := p.adjust(pval, method = "BH"),
           by = .(method, fdr_family)]

  # Keep DICE replication associations significant after FDR correction.
  dice_significant <- dice_all[
    is.finite(FDR) & FDR < 0.05
  ][order(FDR, pval)]

  # Save FDR-significant DICE replication results.
  fwrite(
    dice_significant,
    file.path(
      RESULTS_DICE_DIR,
      "DICE_replication_significant_results.txt"
    ),
    sep = "\t"
  )
}
if (length(log_rows) > 0) {
  fwrite(
    rbindlist(log_rows, fill = TRUE),
    file.path(
      RESULTS_DICE_DIR,
      "DICE_replication_log.txt"
    ),
    sep = "\t"
  )
}

