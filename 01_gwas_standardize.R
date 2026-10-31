# =============================================================================
# GWAS summary-statistics standardization
# =============================================================================
# Purpose:
#   This script standardizes GWAS summary-statistics files for downstream
#   Mendelian randomization and colocalization analyses.
#
# Main steps:
#   1. Read each trait-specific GWAS file from DATA_DIR.
#   2. Detect required columns using a predefined set of aliases.
#   3. Standardize SNP identifiers, alleles, chromosome values, and statistics.
#   4. Remove rows with missing or invalid values.
#   5. Save one standardized <TRAIT>.txt file for each trait.
#
# Input:
#   Raw GWAS files named <TRAIT>.txt under DATA_DIR.
#   The trait list and project directories are loaded from config.R.
#
# Output columns:
#   SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, and pval.
#
# Error handling:
#   If one trait cannot be processed, the error is recorded and the script
#   continues with the remaining traits.
# =============================================================================

source("config.R")
library(data.table)

# Automatically match column names
get_col <- function(dat, aliases, required = TRUE) {
  nm <- tolower(gsub("[^[:alnum:]]", "", names(dat)))
  aliases <- tolower(gsub("[^[:alnum:]]", "", aliases))
  idx <- match(aliases, nm)
  idx <- idx[!is.na(idx)]
  if (!length(idx)) {
    if (required) stop("Column not found: ", paste(aliases, collapse = "/"))
    return(NULL)
  }
  names(dat)[idx[1]]
}
# Record traits that failed processing
failed_traits <- character(0)
for (trait in all_traits) {
  tryCatch({
    message("\nProcessing: ", trait)
    # Automatically locate the raw file
    files <- file.path(DATA_DIR,paste0(trait, c(".txt", ".gz")))
    files <- files[file.exists(files)]
    if (!length(files)) {stop("Raw file not found")}
    gwas <- fread(files[1])
    # Identify basic columns
    snp <- get_col(gwas, c("SNP", "rsids", "variant_id"))
    chr <- get_col(gwas, c("CHR", "chromosome", "chr"))
    pos <- get_col(gwas, c("POS", "GENPOS", "bpos"))
    ea  <- get_col(gwas, c("effect_allele", "A1", "ALLELE1"))
    oa  <- get_col(gwas, c("other_allele",  "A2", "ALLELE2"))
    # Identify beta, OR, SE, and P-value columns
    beta_col <- get_col(gwas, c("effect", "beta", "mtag_beta"), required = FALSE)
    or_col   <- get_col(gwas, c("OR", "odds_ratio"), required = FALSE)
    se_col   <- get_col(gwas, c("SE", "sebeta", "mtag_se"), required = FALSE)
    eaf_col  <- get_col(gwas, c("EAF", "A1FREQ", "freq"), required = FALSE)
    p_col    <- get_col(gwas, c("P", "PVAL", "pval"), required = FALSE)
    # If no standard P-value column exists, look for a log10P column
    if (!is.null(p_col)) {
      pval <- suppressWarnings(as.numeric(gwas[[p_col]]))
      message("Detected P-value column: ", p_col, "; no conversion applied")
    } else {
      log10p_col <- get_col(gwas,c("neglog10p", "minuslog10p", "-log10P","log10P", "log10_p", "mlogp","LOG10P","mtag_pval"),required = FALSE)
      if (!is.null(log10p_col)) {
        log10p <- suppressWarnings(as.numeric(gwas[[log10p_col]]))
        # Treat the input as -log10(P) by default
        pval <- 10^(-log10p)
        message("Detected log10P column: ", log10p_col, "; converted from -log10(P) to P")
      } else {
        stop("No P-value or log10P column found")
      }
    }
    # 4. Get beta; convert log(OR) when only OR is available
    if (!is.null(beta_col)) {
      beta <- suppressWarnings(as.numeric(gwas[[beta_col]]))
    } else if (!is.null(or_col)) {
      OR <- suppressWarnings(as.numeric(gwas[[or_col]]))
      beta <- ifelse(is.finite(OR) & OR > 0, log(OR), NA_real_)
      message("Detected OR; converted to beta = log(OR)")
    } else {
      stop("No beta or OR column found")
    }
    # 5. Get SE; infer it from beta and the P value when unavailable
    if (!is.null(se_col)) {
      se <- suppressWarnings(as.numeric(gwas[[se_col]]))
    } else {
      se <- rep(NA_real_, length(beta))
      valid <- is.finite(beta) & beta != 0 &
        is.finite(pval) & pval > 0 & pval < 1
      z_abs <- qnorm(pval[valid] / 2, lower.tail = FALSE)
      valid_idx <- which(valid)
      good <- is.finite(z_abs) & z_abs > 0
      se[valid_idx[good]] <- abs(beta[valid_idx[good]]) / z_abs[good]
      message("SE column not found; inferred SE from beta and P values")
    }
    
    # 6. Standardize the data format
    out <- data.table(
      SNP = as.character(gwas[[snp]]),
      chr = gwas[[chr]],
      pos = suppressWarnings(as.numeric(gwas[[pos]])),
      effect_allele = toupper(as.character(gwas[[ea]])),
      other_allele = toupper(as.character(gwas[[oa]])),
      eaf = if (!is.null(eaf_col)) {
        suppressWarnings(as.numeric(gwas[[eaf_col]]))
      } else {
        rep(NA_real_, nrow(gwas))
      },
      beta = beta,
      se = se,
      pval = pval
    )
    # 7. Clean SNP identifiers and chromosome values
    out[, SNP := sub(",.*$", "", SNP)]
    out[, SNP := sub("^RS", "rs", SNP, ignore.case = TRUE)]
    out[, chr := sub("^chr", "", trimws(as.character(chr)),
                     ignore.case = TRUE)]
    out <- out[
      !is.na(SNP) & SNP != "" & SNP != "." &
        !is.na(chr) & is.finite(pos) &
        !is.na(effect_allele) & effect_allele != "" &
        !is.na(other_allele) & other_allele != "" &
        is.finite(beta) &
        is.finite(se) & se > 0 &
        is.finite(pval) & pval >= 0 & pval <= 1
    ]
    # 8. Save each trait
    fwrite(out,file.path(DATA_DIR, paste0(trait, ".txt")),sep = "\t")
  }, error = function(e) {
    # Record the failed trait and the error message
    failed_traits <<- c(failed_traits,paste0(trait, ": ", conditionMessage(e)))
    message("[error] Failed to process ", trait, "; skipped. Reason: ", conditionMessage(e))
  })
  # Continue automatically with the next trait after tryCatch
}

# Print the final failure summary after all loops finish
cat("\n\n========== Processing complete ==========\n")

if (length(failed_traits) > 0) {
  cat("The following traits failed (", length(failed_traits), "):\n", sep = "")
  cat(paste0(" - ", failed_traits, collapse = "\n"), "\n")
} else {
  cat("All traits processed successfully!\n")
}
