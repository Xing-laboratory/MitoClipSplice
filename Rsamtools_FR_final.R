# ============================================================================
# Script: Rsamtools_FR_github.R
# Description:
#   This script processes paired-end BAM files (aligned to mitochondrial genome)
#   to extract soft-clip information, compute quality metrics, and generate
#   two tab-delimited output files per input BAM:
#     - Hwenjian: for reads with flags 99 (R1) and 147 (R2)
#     - Lwenjian: for reads with flags 163 (R2) and 83 (R1)
#   The outputs contain soft-clip positions, sequences, base qualities, GC
#   content, and fragment start/end positions.
# Input:
#   "_finally.mt.bam" and "_finally.mt.bam.bai"
# Output:
#   "_Hwenjian.txt" and "_Lwenjian.txt"
# Usage: Rscript Rsamtools_FR_github.R
# ============================================================================




setwd("D:/run/mtbam_mtbambai/")
library(Rsamtools)
library(stringr)

