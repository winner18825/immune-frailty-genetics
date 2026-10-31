# =============================================================================
# config.R -- Central project configuration
# =============================================================================
# Edit the paths and parameters below once to match your environment.
# Each analysis script should begin with source("config.R") to load these settings.
# Run all scripts from the repository root directory.
# =============================================================================

## ---- Project directories ----------------------------------------------------
# Directory containing the harmonized GWAS summary statistics.
# Expected file-naming conventions:
#   <TRAIT>.txt      GWAS summary statistics for MR and colocalization
#   <TRAIT>.txt.gz   eQTL data from Onek1k
#   <TRAIT>.vcf.gz   eQTL data from DICE
DATA_DIR    <- "data/gwas"      # GWAS summary statistics and eQTL data (Onek1k and DICE)
REF_DIR     <- "reference"      # LD reference panel files
RESULTS_DIR <- "/results"       # Directory for all analysis results

## ---- Trait definitions ------------------------------------------------------
# Traits analyzed as outcomes:
all_traits <- c("Depressive_symptoms","Exhaustion","Low_Grip_strength",
                        "Low_Physical_activity","Reaction_time","Slow_Walking_speed",
                        "Social_interaction_and_isolation","Trouble_falling_asleep","Weight_loss")
# Immune cell types used as exposure traits.
cell <- c("cd4et", "cd4nc", "cd4sox4", "cd8et", "cd8nc", "cd8s100b", "monoc", 
                "dc", "bmem", "bin", "nk", "nkr", "mononc", "plasma")
DICE<- c("CD4_NAIVE", "CD8_NAIVE", "MONOCYTES", "B_CELL_NAIVE", "NK", "M2")
# ---- Expected GWAS column names ---------------------------------------------
# The <TRAIT>.txt files should contain the following columns.
# Update these values if your files use different column names.
COL <- list(
  SNP = "SNP", CHR = "chr", BP = "pos",
  EA  = "effect_allele", OA = "other_allele",
  EAF = "eaf", BETA = "beta", SE = "se", P = "pval",
  N="N"
)

# ---- External command-line tools --------------------------------------------
# Path to the PLINK 1.9 executable: https://www.cog-genomics.org/plink/
PLINK_BIN <- "tools/plink/plink"
