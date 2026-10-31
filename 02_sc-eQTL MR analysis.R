# =============================================================================
# Two-sample Mendelian randomization
# =============================================================================
# Exposure: immune-cell eQTL gene expression
# Outcome: multidimensional frailty GWAS traits
#
# The script runs in three stages:
#   1. Select genome-wide significant eQTL SNPs, LD-clump them, and retain
#      instruments with F > 10.
#   2. Run MR, heterogeneity, and MR-Egger intercept tests for one
#      cell/gene/outcome combination.
#   3. Iterate over every available instrument file and outcome trait, then
#      write both per-pair and collated results.
# =============================================================================
source("config.R")
library(TwoSampleMR)
library(ieugwasr)
library(data.table)

## ---- Parameters -------------------------------------------------------------
P_TH       <- 5e-8
CLUMP_KB   <- 10000
CLUMP_R2   <- 0.001
F_FILTER   <- 10
N_DIST     <- 1000

LD_BFILE   <- file.path(REF_DIR, "1000G_EUR", "EUR")
mr_dir     <- file.path(RESULTS_DIR, "mr")
inst_dir   <- file.path(mr_dir, "instruments")
dir.create(mr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(inst_dir, recursive = TRUE, showWarnings = FALSE)
ld_files <- paste0(LD_BFILE, c(".bed", ".bim", ".fam"))

## ---- Helpers ----------------------------------------------------------------
eqtl_required <- c("GENE", "RSID", "beta", "se", "A2", "A1",
                   "A2_FREQ_ONEK1K", "P_VALUE", "N", "CHR", "POS")

## ---- Local LD clumping ------------------------------------------------------
clump_local <- function(dat) {
  clump_input <- data.frame(
    rsid = dat$SNP,
    pval = dat$pval.exposure,
    id = dat$id.exposure,
    stringsAsFactors = FALSE
  )
  pruned <- ieugwasr::ld_clump(
    dat = clump_input,
    clump_kb = CLUMP_KB,
    clump_r2 = CLUMP_R2,
    clump_p = 1,
    bfile = LD_BFILE,
    plink_bin = PLINK_BIN
  )
  dat[dat$SNP %in% pruned$rsid, , drop = FALSE]
}

# =============================================================================
# PART 1 | Instrument preparation
# =============================================================================
message("\n========== Part 1: instrument preparation ==========\n")

for (cellid in cell) {
  message(">>> exposure: ", cellid)
  input_path <- file.path(DATA_DIR, paste0(cellid, "_eqtl_table.txt.gz"))
  if (!file.exists(input_path)) {
    message("  [skip] file not found: ", input_path)
    next
  }
  
  exposure_raw <- tryCatch(
    fread(input_path, sep = "\t"),
    error = function(e) {
      message("  [error] reading file: ", e$message)
      NULL
    }
  )
  if (is.null(exposure_raw)) next

  inst_dir_cell <- file.path(inst_dir, cellid)
  dir.create(inst_dir_cell, recursive = TRUE, showWarnings = FALSE)
  genes <- unique(as.character(exposure_raw$GENE))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  
  for (geneid in genes) {
    gene_raw <- exposure_raw[GENE == geneid]
    exposure_dat <- tryCatch(
      format_data(
        gene_raw,
        type = "exposure",
        snp_col = "RSID",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "A2",
        other_allele_col = "A1",
        eaf_col = "A2_FREQ_ONEK1K",
        pval_col = "P_VALUE",
        samplesize_col = "N",
        chr_col = "CHR",
        pos_col = "POS",
        phenotype_col = "GENE"
      ),
      error = function(e) {
        message("  [error] ", geneid, ": ", e$message)
        NULL
      }
    )
    if (is.null(exposure_dat)) next
    exposure_dat$exposure <- paste(cellid, geneid, sep = ":")
    exposure_dat$id.exposure <- paste(cellid, geneid, sep = ":")
    
    sig <- exposure_dat[is.finite(exposure_dat$pval.exposure) &
                          exposure_dat$pval.exposure < P_TH, , drop = FALSE]
    sig <- sig[!duplicated(sig$SNP), , drop = FALSE]
    if (nrow(sig) == 0) {
      message("  [skip] ", geneid, ": no significant SNP")
      next
    }
    message("  ", geneid, ": significant SNPs = ", nrow(sig))
    
    clumped <- tryCatch(
      clump_local(sig),
      error = function(e) {
        message("  [error] clumping ", geneid, ": ", e$message)
        NULL
      }
    )
    if (is.null(clumped) || nrow(clumped) == 0) {
      message("  [skip] ", geneid, ": no SNP after clumping")
      next
    }
    message("  ", geneid, ": after clumping = ", nrow(clumped))
    
    clumped$R2 <- with(clumped, 2 * beta.exposure^2 *
                         eaf.exposure * (1 - eaf.exposure))
    clumped$F <- with(clumped, (samplesize.exposure - 2) * R2 / (1 - R2))
    out_tab <- clumped[is.finite(clumped$F) & clumped$F > F_FILTER, , drop = FALSE]
    if (nrow(out_tab) == 0) {
      message("  [skip] ", geneid, ": no SNP after F filter")
      next
    }
    message("  ", geneid, ": after F filter = ", nrow(out_tab))
    
    fwrite(out_tab,
           file.path(inst_dir_cell, paste0(geneid, ".instruments.txt")),
           sep = "\t")
  }
}

# =============================================================================
# PART 2 | MR analysis for one cell -> gene -> outcome combination
# =============================================================================
run_one_MR <- function(cellid, geneid, outcome_trait) {
  message("\n------ MR: ", cellid, ":", geneid, " -> ", outcome_trait, " ------")
  exposure_file <- file.path(inst_dir, cellid,
                             paste0(geneid, ".instruments.txt"))
  outcome_file <- file.path(DATA_DIR, paste0(outcome_trait, ".txt"))
  if (!file.exists(exposure_file) || !file.exists(outcome_file)) {
    message("  [skip] input file not found")
    return(NULL)
  }
  
  pair_dir <- file.path(mr_dir, cellid, geneid)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  
  exposure_dat <- tryCatch(
    read_exposure_data(
      filename = exposure_file, sep = "\t",
      snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure", pval_col = "pval.exposure",
      samplesize_col = "samplesize.exposure",
      chr_col = "chr.exposure", pos_col = "pos.exposure", clump = FALSE
    ),
    error = function(e) {
      message("  [error] exposure: ", e$message)
      NULL
    }
  )
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) return(NULL)
  
  outcome_dat <- tryCatch(
    read_outcome_data(
      snps = exposure_dat$SNP, filename = outcome_file, sep = "\t",
      snp_col = COL$SNP, beta_col = COL$BETA, se_col = COL$SE,
      effect_allele_col = COL$EA, other_allele_col = COL$OA,
      eaf_col = COL$EAF, pval_col = COL$P,
      chr_col = COL$CHR, pos_col = COL$BP
    ),
    error = function(e) {
      message("  [error] outcome: ", e$message)
      NULL
    }
  )
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    message("  [skip] no exposure SNPs found in outcome")
    return(NULL)
  }
  
  exposure_dat$exposure <- exposure_dat$id.exposure <- paste(cellid, geneid, sep = ":")
  outcome_dat$outcome <- outcome_dat$id.outcome <- outcome_trait
  dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
  dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
  if (nrow(dat) == 0) {
    message("  [skip] no SNPs after harmonisation")
    return(NULL)
  }
  message("  SNPs after harmonisation: ", nrow(dat))
  
  mr_res <- tryCatch(mr(dat), error = function(e) {
    message("  [error] MR: ", e$message)
    NULL
  })
  if (!is.null(mr_res) && nrow(mr_res) > 0) {
    mr_res$cell <- cellid
    mr_res$gene <- geneid
    mr_res$outcome_trait <- outcome_trait
    mr_res$OR <- exp(mr_res$b)
    mr_res$OR_lci95 <- exp(mr_res$b - 1.96 * mr_res$se)
    mr_res$OR_uci95 <- exp(mr_res$b + 1.96 * mr_res$se)
    mr_res$beta_lci95 <- mr_res$b - 1.96 * mr_res$se
    mr_res$beta_uci95 <- mr_res$b + 1.96 * mr_res$se
    fwrite(mr_res, file.path(pair_dir, "MR_results.txt"), sep = "\t")
  }
  
  het_res <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(het_res) && nrow(het_res) > 0) {
    het_res$cell <- cellid
    het_res$gene <- geneid
    het_res$outcome_trait <- outcome_trait
    fwrite(het_res, file.path(pair_dir, "heterogeneity.txt"), sep = "\t")
  }
  
  plt_res <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(plt_res) && nrow(plt_res) > 0) {
    plt_res$cell <- cellid
    plt_res$gene <- geneid
    plt_res$outcome_trait <- outcome_trait
    fwrite(plt_res, file.path(pair_dir, "pleiotropy.txt"), sep = "\t")
  }
  
  presso_result <- NULL
  if (nrow(dat) >= 4 && requireNamespace("MRPRESSO", quietly = TRUE)) {
    presso_result <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = as.data.frame(dat), NbDistribution = N_DIST,
        SignifThreshold = 0.05
      ),
      error = function(e) {
        message("  [warning] MR-PRESSO: ", e$message)
        NULL
      }
    )
    if (!is.null(presso_result)) {
      saveRDS(presso_result, file.path(pair_dir, "MR-PRESSO.rds"))
    }
  }
  
  message("  [done] ", cellid, ":", geneid, " -> ", outcome_trait)
  list(mr_res = mr_res, het_res = het_res, plt_res = plt_res,
       presso = presso_result, nsnp = nrow(dat), cell = cellid,
       gene = geneid, outcome_trait = outcome_trait)
}

# =============================================================================
# PART 3 | Run every available pair and collate results
# =============================================================================
message("\n========== Part 3: MR for all pairs ==========\n")

all_mr <- list()
all_het <- list()
all_plt <- list()
log_rows <- list()

for (cellid in cell) {
  instrument_files <- list.files(
    file.path(inst_dir, cellid),
    pattern = "\\.instruments\\.txt$",
    full.names = FALSE
  )
  if (length(instrument_files) == 0) {
    message(">>> ", cellid, ": no instrument files")
    next
  }
  for (instrument_file in instrument_files) {
    geneid <- sub("\\.instruments\\.txt$", "", instrument_file)
    for (outcome_trait in all_traits) {
      result <- run_one_MR(cellid, geneid, outcome_trait)
      if (is.null(result)) next
      if (!is.null(result$mr_res) && nrow(result$mr_res) > 0)
        all_mr[[length(all_mr) + 1L]] <- result$mr_res
      if (!is.null(result$het_res) && nrow(result$het_res) > 0)
        all_het[[length(all_het) + 1L]] <- result$het_res
      if (!is.null(result$plt_res) && nrow(result$plt_res) > 0)
        all_plt[[length(all_plt) + 1L]] <- result$plt_res
      log_rows[[length(log_rows) + 1L]] <- data.table(
        cell = cellid, gene = geneid, outcome_trait = outcome_trait,
        nsnp = result$nsnp, status = "ok"
      )
    }
  }
}

if (length(log_rows) > 0)
  fwrite(rbindlist(log_rows, fill = TRUE),
         file.path(mr_dir, "MR_analysis_log.txt"), sep = "\t")

if (length(all_mr) > 0) {
  df_all <- rbindlist(all_mr, fill = TRUE)
  # Save all MR results for the selected methods.
  df_mr <- df_all[
    method %in% c("Inverse variance weighted", "Wald ratio") &
      is.finite(pval)
  ]
  fwrite(df_mr, file.path(mr_dir, "MR_all_methods_results.txt"),sep = "\t")
  # Define the five physical frailty components as one testing family.

  physical_frailty_traits <- c(
    "Weight_loss",
    "Exhaustion",
    "Low_Physical_activity",
    "Slow_Walking_speed",
    "Low_Grip_strength"
  )
  # Assign the FDR testing family.
  # The five physical frailty components share one family.
  # Each remaining outcome is corrected separately.
  df_mr[, fdr_family := fifelse(
    outcome_trait %in% physical_frailty_traits,
    "physical_frailty",
    outcome_trait
  )]

  # Apply the Benjamini-Hochberg correction separately for each MR method
  # and each predefined testing family.
  df_mr[, FDR := p.adjust(pval, method = "BH"),
        by = .(method, fdr_family)]

  # Keep associations significant after FDR correction.
  df_significant <- df_mr[
    is.finite(FDR) & FDR < 0.05
  ][order(FDR, pval)]

  # Save the FDR-significant MR results for colocalization analysis.
  fwrite(df_significant, file.path(mr_dir, "MR_significant_results.txt"),sep = "\t")

}
if (length(all_het) > 0)
  fwrite(rbindlist(all_het, fill = TRUE),
         file.path(mr_dir, "MR_heterogeneity_results.txt"), sep = "\t")
if (length(all_plt) > 0)
  fwrite(rbindlist(all_plt, fill = TRUE),
         file.path(mr_dir, "MR_pleiotropy_results.txt"), sep = "\t")

cat("\nMR complete. Results in: ", normalizePath(mr_dir, mustWork = FALSE), "\n", sep = "")
