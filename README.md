# Immune Cell eQTL Mendelian Randomization, Colocalization, and DICE Replication

Analysis code for assessing the causal effects of immune-cell-specific gene
expression on multidimensional frailty traits, followed by Bayesian
colocalization and replication using DICE immune-cell eQTL data.

The workflow uses publicly available GWAS summary statistics, OneK1K
single-cell eQTL data, DICE eQTL data, local LD reference data, and open-source
R packages.

---

## Analysis pipeline

The analysis is organized into four sequential modules. Each script reads the
shared configuration from `config.R` and writes results to `RESULTS_DIR`.

| Step | Script | Analysis | Main output |
|------|--------|----------|-------------|
| 1 | `01_gwas_standardize.R` | Standardize GWAS summary statistics | `data/gwas/<TRAIT>.txt` |
| 2 | `02_sc-eQTL MR analysis.R` | Two-sample MR using OneK1K sc-eQTL instruments | `results/mr/` |
| 3 | `03_Colocalization analysis.R` | Bayesian colocalization of MR-positive pairs with all frailty traits | `results/coloc/` |
| 4 | `04_DICE_replication_MR.R` | DICE-based replication MR | `results/dice_replication/` |

The overall workflow is:

```text
Raw GWAS data
    -> GWAS standardization
    -> OneK1K eQTL instrument preparation
    -> Two-sample MR
    -> MR-positive cell-gene pairs
    -> Colocalization with frailty GWAS traits
    -> DICE replication MR
```

---

## Software requirements

Use R (recommended version >= 4.1) with the following packages:

```r
install.packages(c(
  "data.table",
  "TwoSampleMR",
  "ieugwasr",
  "coloc"
))
```

Optional packages:

```r
install.packages(c("MRPRESSO", "locuscomparer"))
```

`MRPRESSO` is used only when available in the MR script. `locuscomparer` can
be used for optional colocalization plots.

The MR instrument-preparation step also requires:

- PLINK 1.9
- A 1000 Genomes EUR LD reference panel containing `.bed`, `.bim`, and `.fam`
  files

Set the PLINK executable and LD reference locations in `config.R`.

---

## Configuration

Edit `config.R` before running the analysis. The main settings are:

- `DATA_DIR`: directory containing GWAS, OneK1K eQTL, and DICE input files
- `REF_DIR`: directory containing the LD reference panel
- `RESULTS_DIR`: directory where all analysis results are written
- `all_traits`: the nine frailty-related outcome traits
- `cell`: the OneK1K immune-cell labels
- `DICE`: the DICE immune-cell labels
- `COL`: standardized GWAS column names
- `PLINK_BIN`: path to the PLINK 1.9 executable

The current DICE mapping is:

| OneK1K / colocalization label | DICE label |
|-------------------------------|------------|
| `cd4nc` | `CD4_NAIVE` |
| `cd8nc` | `CD8_NAIVE` |
| `monoc` | `MONOCYTES` |
| `bin` | `B_CELL_NAIVE` |
| `nk` | `NK` |
| `mononc` | `M2` |

The cell labels in colocalization results must remain the OneK1K labels. The
DICE replication script applies the mapping above when locating DICE files.

---

## Input data

GWAS summary statistics are not included in this repository. Place the raw
files under `DATA_DIR` using the trait names defined in `config.R`.

### Raw GWAS files

The standardization script searches for:

```text
<TRAIT>.txt
<TRAIT>.gz
```

After standardization, the outcome files should contain:

```text
SNP, chr, pos, effect_allele, other_allele,
eaf, beta, se, pval, N
```

The `N` column is required by the colocalization script. If the current GWAS
standardization script does not carry `N` into the output, retain or add this
column before running colocalization.

### OneK1K sc-eQTL files

Store one file per immune cell under `DATA_DIR`:

```text
<cell>_eqtl_table.txt.gz
```

Required columns:

```text
GENE, RSID, beta, se, A2, A1,
A2_FREQ_ONEK1K, P_VALUE, N, CHR, POS
```

### DICE eQTL files

Store one compressed VCF-style file per DICE cell under `DATA_DIR`:

```text
CD4_NAIVE.vcf.gz
CD8_NAIVE.vcf.gz
MONOCYTES.vcf.gz
B_CELL_NAIVE.vcf.gz
NK.vcf.gz
M2.vcf.gz
```

Required DICE columns:

```text
#CHROM, POS, ID, REF, ALT, Gene,
GeneSymbol, Pvalue, Beta, Cell
```

The DICE replication script derives the standard error from `Beta` and the
two-sided `Pvalue`, because the DICE input format does not provide an SE
column.

### LD reference files

The MR script expects a PLINK fileset with this prefix:

```text
<REF_DIR>/1000G_EUR/EUR.bed
<REF_DIR>/1000G_EUR/EUR.bim
<REF_DIR>/1000G_EUR/EUR.fam
```

---

## How to run

Run all scripts from the repository root so that `source("config.R")` resolves
correctly.

### 1. Standardize GWAS summary statistics

```r
source("01_gwas_standardize.R")
```

This creates standardized `<TRAIT>.txt` files under `DATA_DIR`.

### 2. Run OneK1K sc-eQTL MR

```r
source("02_sc-eQTL MR analysis.R")
```

The script selects genome-wide significant eQTL variants using `P < 5e-8`,
performs local LD clumping, applies an instrument-strength filter of `F > 10`,
and runs MR across all available cell-gene and frailty-trait pairs.

Important output files include:

```text
results/mr/MR_all_methods_results.txt
results/mr/MR_IVW_results_summary.txt
results/mr/MR_IVW_significant_results.txt
results/mr/MR_heterogeneity_results.txt
results/mr/MR_pleiotropy_results.txt
```

### 3. Run Bayesian colocalization

```r
source("03_Colocalization analysis.R")
```

The script reads the MR-positive results, keeps significant inverse-variance-
weighted and Wald-ratio associations, and tests the corresponding cell-gene
pairs against all nine traits.

The outcome models are:

- `quant`: `Depressive_symptoms` and `Reaction_time`
- `cc`: `Exhaustion`, `Low_Grip_strength`, `Low_Physical_activity`,
  `Slow_Walking_speed`, `Social_interaction_and_isolation`,
  `Trouble_falling_asleep`, and `Weight_loss`

For case-control traits, the script uses the trait-specific case fractions
specified in `CASE_FRACTION`.

Main output files include:

```text
results/coloc/coloc_summary_all.tsv
results/coloc/<cell>__<gene>__<trait>.summary.tsv
results/coloc/<cell>__<gene>__<trait>.snp.tsv
results/coloc/coloc_errors.tsv
```

### 4. Run DICE replication MR

```r
source("04_DICE_replication_MR.R")
```

The script selects colocalization-positive pairs with `PP.H4 >= 0.90`, maps
the OneK1K cell labels to DICE labels, selects DICE eQTL variants with
`Pvalue < 5e-8`, harmonizes them with the corresponding GWAS, and runs
replication MR.

Main output files include:

```text
results/dice_replication/DICE_replication_all_results.txt
results/dice_replication/DICE_replication_log.txt
results/dice_replication/<cell>/<gene>/<trait>.MR.txt
```

---

## Reproducibility notes

- Use the same trait names in `config.R`, GWAS filenames, and the colocalization
  results.
- Keep SNP identifiers consistent across OneK1K, DICE, and GWAS files.
- Check effect-allele orientation and allele frequencies before interpreting
  harmonized MR results.
- The DICE replication step uses only the six mapped immune-cell populations
  listed in `DICE_MAP`.
- The scripts write intermediate and final results under `RESULTS_DIR`; make
  sure this directory exists and is writable.

---

## Method references

- **TwoSampleMR** — Hemani et al. *eLife* 2018.
- **MR-PRESSO** — Verbanck et al. *Nature Genetics* 2018.
- **coloc** — Giambartolomei et al. *PLoS Genetics* 2014.
- **PLINK** — Purcell et al. *American Journal of Human Genetics* 2007.

## Citation

If you use this code, please cite the associated study and the software listed
above.

## License

Add the project license here before public release.
