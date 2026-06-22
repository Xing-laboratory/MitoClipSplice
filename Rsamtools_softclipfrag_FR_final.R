# ==============================================================================
# Script: Rsamtools_softclipfrag_FR_final.R
# Description:
#   This script processes the *_Hwenjian.txt and *_Lwenjian.txt files generated
#   by Rsamtools_FR_final.R. It compares fragments that contain
#   soft-clipped reads (H1/L1) versus those without (H0/L0) across three
#   metrics: mapping quality (Q30 proportion), GC content (median), and insert
#   size (median). For each input file, it generates distribution tables and a
#   summary result file.
#   The script is split into two independent sections:
#     - Part A: Process Hwenjian files (proper pairs flags 99/147)
#     - Part B: Process Lwenjian files (proper pairs flags 163/83)
# Input:
#   "_Hwenjian.txt" and "_Lwenjian.txt"
# Output:
#   For each input file, output files are generated:
#     - {prefix}_H1_Mqmean_table.txt / - {prefix}_H0_Mqmean_table.txt
#     - {prefix}_H1_GC_table.txt / - {prefix}_H0_GC_table.txt
#     - {prefix}_H1_insert_table.txt / - {prefix}_H0_insert_table.txt
#     - {prefix}_H_result.txt
#     - {prefix}_L1_Mqmean_table.txt / - {prefix}_L0_Mqmean_table.txt
#     - {prefix}_L1_GC_table.txt / - {prefix}_L0_GC_table.txt
#     - {prefix}_L1_insert_table.txt / - {prefix}_L0_insert_table.txt
#     - {prefix}_L_result.txt
# Usage: Rscript Rsamtools_softclipfrag_FR_final.R
# ==============================================================================




# ------------------------------------------------------------------------------
# Part 1: Process H-strand files (flags 99 & 147)
# ------------------------------------------------------------------------------


# 1.1 Set paths for H-strand files
root_path <- "D:/run/"
input_subdir <- "Hwenjian/"
input_full_path <- paste0(root_path, input_subdir)
output_parent_subdir <- "softclipfrag/"
output_root_path <- paste0(root_path, output_parent_subdir)

# Define output subdirectories for H-strand
output_subdirs <- c(
  "H1_Mqmean_table",
  "H0_Mqmean_table",
  "H1_GC_table",
  "H0_GC_table",
  "H1_insert_table",
  "H0_insert_table",
  "H_result"
)

# 1.2 Create all output directories
for (subdir in output_subdirs) {
  dir_path <- paste0(output_root_path, subdir, "/")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# 1.3 Get list of H-strand files
ff1 <- list.files(input_full_path, pattern = "\\_Hwenjian.txt$")
cat("Number of H-strand files found:", length(ff1), "\n")
print(ff1)

# 1.4 Process each H-strand file
for(i in 1:length(ff1)){
  # Read data
  H <- read.table(paste0(input_full_path, ff1[i]), sep = "\t", header = TRUE)
  num_H <- nrow(H)
  
  # Split into soft-clipped (H1) and non-soft-clipped (H0) groups
  H1 <- subset(H, H$R1S5type != 0 | H$R2S5type != 0)
  number_H1 <- nrow(H1)
  number_H1_pro <- number_H1 / num_H
  
  H0 <- subset(H, H$R1S5type == 0 & H$R2S5type == 0)
  number_H0 <- nrow(H0)
  number_H0_pro <- number_H0 / num_H
  
  # -------------------- Mean mapping quality (Mqmean) analysis --------------------
  H1_Mqmean_table <- as.data.frame(table(
    cut(na.omit(c(H1$R1Mqmean, H1$R2Mqmean)), 
        breaks = seq(0, 50, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  H1_Mqmean_table[,3] <- H1_Mqmean_table[,2] / (number_H1 * 2)
  colnames(H1_Mqmean_table) <- c("Type","Number","Proportion")
  
  H0_Mqmean_table <- as.data.frame(table(
    cut(na.omit(c(H0$R1Mqmean, H0$R2Mqmean)), 
        breaks = seq(0, 50, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  H0_Mqmean_table[,3] <- H0_Mqmean_table[,2] / (number_H0 * 2)
  colnames(H0_Mqmean_table) <- c("Type","Number","Proportion")
  
  # Q30 proportion
  H1_Mqmean_pro <- sum(H1_Mqmean_table[31:50,3]) / sum(H1_Mqmean_table[,3])
  H0_Mqmean_pro <- sum(H0_Mqmean_table[31:50,3]) / sum(H0_Mqmean_table[,3])
  
  # Write Mqmean tables
  write.table(
    x = H1_Mqmean_table, 
    file = paste0(output_root_path, "H1_Mqmean_table/", sub("_Hwenjian.txt", "_H1_Mqmean_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = H0_Mqmean_table, 
    file = paste0(output_root_path, "H0_Mqmean_table/", sub("_Hwenjian.txt", "_H0_Mqmean_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- GC content analysis --------------------
  H1_GC_table <- as.data.frame(table(
    cut(na.omit(c(H1$R1GC, H1$R2GC)), 
        breaks = seq(0, 100, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  H1_GC_table[,3] <- H1_GC_table[,2] / (number_H1 * 2)
  colnames(H1_GC_table) <- c("Type","Number","Proportion")
  
  H0_GC_table <- as.data.frame(table(
    cut(na.omit(c(H0$R1GC, H0$R2GC)), 
        breaks = seq(0, 100, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  H0_GC_table[,3] <- H0_GC_table[,2] / (number_H0 * 2)
  colnames(H0_GC_table) <- c("Type","Number","Proportion")
  
  # Median GC
  H1_GC_median <- median(c(H1$R1GC, H1$R2GC), na.rm = TRUE)
  H0_GC_median <- median(c(H0$R1GC, H0$R2GC), na.rm = TRUE)
  
  # Write GC tables
  write.table(
    x = H1_GC_table, 
    file = paste0(output_root_path, "H1_GC_table/", sub("_Hwenjian.txt", "_H1_GC_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = H0_GC_table, 
    file = paste0(output_root_path, "H0_GC_table/", sub("_Hwenjian.txt", "_H0_GC_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- Insert size analysis --------------------
  H1_insert_table <- as.data.frame(table(
    cut(na.omit(H1$R1isize), 
        breaks = seq(0, 600, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  H1_insert_table[,3] <- H1_insert_table[,2] / number_H1
  colnames(H1_insert_table) <- c("Type","Number","Proportion")
  
  H0_insert_table <- as.data.frame(table(
    cut(na.omit(H0$R1isize), 
        breaks = seq(0, 600, by = 1), 
        right = FALSE, include.lowest = TRUE)
  ))
  H0_insert_table[,3] <- H0_insert_table[,2] / number_H0
  colnames(H0_insert_table) <- c("Type","Number","Proportion")
  
  # Median insert size
  H1_insert_median <- median(H1$R1isize, na.rm = TRUE)
  H0_insert_median <- median(H0$R1isize, na.rm = TRUE)
  
  # Write insert size tables
  write.table(
    x = H1_insert_table, 
    file = paste0(output_root_path, "H1_insert_table/", sub("_Hwenjian.txt", "_H1_insert_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = H0_insert_table, 
    file = paste0(output_root_path, "H0_insert_table/", sub("_Hwenjian.txt", "_H0_insert_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- Summary result table --------------------
  H_result <- rbind(
    number_H1_pro, number_H0_pro,
    H1_Mqmean_pro, H0_Mqmean_pro,
    H1_GC_median, H0_GC_median,
    H1_insert_median, H0_insert_median
  )
  H_result <- as.data.frame(H_result)
  rownames(H_result) <- c(
    "Soft-clipped fragment proportion",
    "Non-soft-clipped fragment proportion",
    "Soft-clipped fragment Q30 proportion",
    "Non-soft-clipped fragment Q30 proportion",
    "Soft-clipped fragment GC median",
    "Non-soft-clipped fragment GC median",
    "Soft-clipped fragment insert median",
    "Non-soft-clipped fragment insert median"
  )
  colnames(H_result) <- c("H strand")
  
  # Write summary
  write.table(
    x = H_result, 
    file = paste0(output_root_path, "H_result/", sub("_Hwenjian.txt", "_H_result.txt", ff1[i])), 
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
  )
  
  cat(paste("H-strand file processed:", i, "/", length(ff1), " -", ff1[i], "\n"))
}

cat("All H-strand files processed! Output directory:", output_root_path, "\n")




# ------------------------------------------------------------------------------
# Part 2: Process L-strand files (flags 163 & 83)
# ------------------------------------------------------------------------------


# 2.1 Set paths for L-strand files
root_path <- "D:/run/"
input_subdir <- "Lwenjian/"
input_full_path <- paste0(root_path, input_subdir)
output_parent_subdir <- "softclipfrag/"
output_root_path <- paste0(root_path, output_parent_subdir)

# Define output subdirectories for L-strand
output_subdirs <- c(
  "L1_Mqmean_table",
  "L0_Mqmean_table",
  "L1_GC_table",
  "L0_GC_table",
  "L1_insert_table",
  "L0_insert_table",
  "L_result"
)

# 2.2 Create all output directories
for (subdir in output_subdirs) {
  dir_path <- paste0(output_root_path, subdir, "/")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# 2.3 Get list of L-strand files
ff1 <- list.files(input_full_path, pattern = "\\_Lwenjian.txt$")
cat("Number of L-strand files found:", length(ff1), "\n")
print(ff1)

# 2.4 Process each L-strand file
for(i in 1:length(ff1)){
  # Read data
  L <- read.table(paste0(input_full_path, ff1[i]), sep = "\t", header = TRUE)
  num_L <- nrow(L)
  
  # Split into soft-clipped (L1) and non-soft-clipped (L0) groups
  L1 <- subset(L, L$R1S5type != 0 | L$R2S5type != 0)
  number_L1 <- nrow(L1)
  number_L1_pro <- number_L1 / num_L
  
  L0 <- subset(L, L$R1S5type == 0 & L$R2S5type == 0)
  number_L0 <- nrow(L0)
  number_L0_pro <- number_L0 / num_L
  
  # -------------------- Mean mapping quality (Mqmean) analysis --------------------
  L1_Mqmean_table <- as.data.frame(table(
    cut(na.omit(c(L1$R1Mqmean, L1$R2Mqmean)), 
        breaks = seq(0, 50, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  L1_Mqmean_table[,3] <- L1_Mqmean_table[,2] / (number_L1 * 2)
  colnames(L1_Mqmean_table) <- c("Type","Number","Proportion")
  
  L0_Mqmean_table <- as.data.frame(table(
    cut(na.omit(c(L0$R1Mqmean, L0$R2Mqmean)), 
        breaks = seq(0, 50, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  L0_Mqmean_table[,3] <- L0_Mqmean_table[,2] / (number_L0 * 2)
  colnames(L0_Mqmean_table) <- c("Type","Number","Proportion")
  
  # Q30 proportion
  L1_Mqmean_pro <- sum(L1_Mqmean_table[31:50,3]) / sum(L1_Mqmean_table[,3])
  L0_Mqmean_pro <- sum(L0_Mqmean_table[31:50,3]) / sum(L0_Mqmean_table[,3])
  
  # Write Mqmean tables
  write.table(
    x = L1_Mqmean_table, 
    file = paste0(output_root_path, "L1_Mqmean_table/", sub("_Lwenjian.txt", "_L1_Mqmean_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = L0_Mqmean_table, 
    file = paste0(output_root_path, "L0_Mqmean_table/", sub("_Lwenjian.txt", "_L0_Mqmean_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- GC content analysis --------------------
  L1_GC_table <- as.data.frame(table(
    cut(na.omit(c(L1$R1GC, L1$R2GC)), 
        breaks = seq(0, 100, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  L1_GC_table[,3] <- L1_GC_table[,2] / (number_L1 * 2)
  colnames(L1_GC_table) <- c("Type","Number","Proportion")
  
  L0_GC_table <- as.data.frame(table(
    cut(na.omit(c(L0$R1GC, L0$R2GC)), 
        breaks = seq(0, 100, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  L0_GC_table[,3] <- L0_GC_table[,2] / (number_L0 * 2)
  colnames(L0_GC_table) <- c("Type","Number","Proportion")
  
  # Median GC
  L1_GC_median <- median(c(L1$R1GC, L1$R2GC), na.rm = TRUE)
  L0_GC_median <- median(c(L0$R1GC, L0$R2GC), na.rm = TRUE)
  
  # Write GC tables
  write.table(
    x = L1_GC_table, 
    file = paste0(output_root_path, "L1_GC_table/", sub("_Lwenjian.txt", "_L1_GC_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = L0_GC_table, 
    file = paste0(output_root_path, "L0_GC_table/", sub("_Lwenjian.txt", "_L0_GC_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- Insert size analysis --------------------
  L1_insert_table <- as.data.frame(table(
    cut(na.omit(L1$R1isize), 
        breaks = seq(0, 600, by = 1), 
        right = FALSE, 
        include.lowest = TRUE)
  ))
  L1_insert_table[,3] <- L1_insert_table[,2] / number_L1
  colnames(L1_insert_table) <- c("Type","Number","Proportion")
  
  L0_insert_table <- as.data.frame(table(
    cut(na.omit(L0$R1isize), 
        breaks = seq(0, 600, by = 1), 
        right = FALSE, include.lowest = TRUE)
  ))
  L0_insert_table[,3] <- L0_insert_table[,2] / number_L0
  colnames(L0_insert_table) <- c("Type","Number","Proportion")
  
  # Median insert size
  L1_insert_median <- median(L1$R1isize, na.rm = TRUE)
  L0_insert_median <- median(L0$R1isize, na.rm = TRUE)
  
  # Write insert size tables
  write.table(
    x = L1_insert_table, 
    file = paste0(output_root_path, "L1_insert_table/", sub("_Lwenjian.txt", "_L1_insert_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    x = L0_insert_table, 
    file = paste0(output_root_path, "L0_insert_table/", sub("_Lwenjian.txt", "_L0_insert_table.txt", ff1[i])), 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # -------------------- Summary result table --------------------
  L_result <- rbind(
    number_L1_pro, number_L0_pro,
    L1_Mqmean_pro, L0_Mqmean_pro,
    L1_GC_median, L0_GC_median,
    L1_insert_median, L0_insert_median
  )
  L_result <- as.data.frame(L_result)
  rownames(L_result) <- c(
    "Soft-clipped fragment proportion",
    "Non-soft-clipped fragment proportion",
    "Soft-clipped fragment Q30 proportion",
    "Non-soft-clipped fragment Q30 proportion",
    "Soft-clipped fragment GC median",
    "Non-soft-clipped fragment GC median",
    "Soft-clipped fragment insert median",
    "Non-soft-clipped fragment insert median"
  )
  colnames(L_result) <- c("L strand")
  
  # Write summary
  write.table(
    x = L_result, 
    file = paste0(output_root_path, "L_result/", sub("_Lwenjian.txt", "_L_result.txt", ff1[i])), 
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
  )
  
  cat(paste("L-strand file processed:", i, "/", length(ff1), " -", ff1[i], "\n"))
}

cat("All L-strand files processed! Output directory:", output_root_path, "\n")




