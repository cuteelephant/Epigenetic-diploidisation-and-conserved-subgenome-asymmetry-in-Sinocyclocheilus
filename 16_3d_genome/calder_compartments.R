#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: calder_compartments.R <input.cool> <binsize> <outdir>")
}

input_cool <- args[1]
bin_size   <- as.integer(args[2])
out_dir    <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(CALDER)
})

# CALDER expects per-chromosome dump files in a "raw" format (chr_i, chr_j,
# count). The user can convert .cool with `cooler dump --chrom X --chrom2 X`
# in a wrapper for each chromosome. Here we delegate the work to CALDER's
# convenience function CALDER_main(), which reads the .cool directly.

# Sinocyclocheilus tianlinensis / longibarbatus chromosome list — replace as
# appropriate; pass an env var SINO_CHROMS="tl1,tl2,..." to override.
chrom_list_env <- Sys.getenv("SINO_CHROMS", unset = NA)
if (!is.na(chrom_list_env) && nchar(chrom_list_env) > 0) {
  chrom_names <- unlist(strsplit(chrom_list_env, ","))
} else {
  chrom_names <- paste0("tl", 1:50)
}

CALDER_main(
  contact_file_dump      = input_cool,
  chrs                   = chrom_names,
  bin_size               = bin_size,
  genome                 = "other",
  save_dir               = out_dir,
  save_intermediate_data = FALSE,
  n_cores                = 8,
  sub_domains            = TRUE
)

cat(sprintf("[CALDER] Compartments written to %s\n", out_dir))
