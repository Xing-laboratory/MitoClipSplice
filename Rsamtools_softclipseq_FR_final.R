# 定义输入和输出路径
input_dir <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Hwenjian/"
output_base <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/softclipseq/"

# 创建输出目录（在循环外创建，避免重复创建）
for(dir_name in c("H_Snumber_table", "H_Sbase_table", "H_Sb_table", "H_type3_table")) {
  dir.create(file.path(output_base, dir_name), recursive = TRUE, showWarnings = FALSE)
}

# 获取文件列表
ff1 <- list.files(input_dir, pattern = "_Hwenjian.txt$")

# 处理每个文件
for(i in 1:length(ff1)) {
  # 读取数据
  H <- read.table(paste(input_dir, ff1[i], sep = ""), sep = "\t", head = TRUE)
  num_H <- nrow(H)
  
  # 1. softclip长度分布统计
  H_Snumber_tab <- as.data.frame(table(cut(na.omit(c(H$R1S5number, H$R2S5number)), 
                                           breaks = seq(0, 30, by = 1), 
                                           right = FALSE, include.lowest = TRUE)))
  H_Snumber_table <- H_Snumber_tab[-1,]
  H_Snumber_table[,3] <- H_Snumber_table[,2] / sum(H_Snumber_table[,2])
  colnames(H_Snumber_table) <- c("Type","Number","Proportion")
  
  # 2. softclip序列碱基组成统计
  H_Sbase <- c(H$R1S5base, H$R2S5base)
  H_Sbase_seq <- paste(H_Sbase, collapse = "")
  H_Sbase_table <- table(strsplit(H_Sbase_seq, "")[[1]])
  H_Sbase_table <- as.data.frame(H_Sbase_table)
  H_Sbase_table <- H_Sbase_table[!grepl("N", H_Sbase_table$Var1), ]
  H_Sbase_table[,3] <- H_Sbase_table[,2]/sum(H_Sbase_table[,2])
  colnames(H_Sbase_table) <- c("base1234","H_number","H_proportion")
  
  # 3. softclip第一个碱基组成统计
  H_R11 <- subset(H, H$R1S5type != 0)
  H_R21 <- subset(H, H$R2S5type != 0)
  H_Sb <- c(H_R11$R1S5b1, H_R21$R2S5b1)
  H_Sb_seq <- paste(H_Sb, collapse = "")
  H_Sb_table <- table(strsplit(H_Sb_seq, "")[[1]])
  H_Sb_table <- as.data.frame(H_Sb_table)
  H_Sb_table <- H_Sb_table[!grepl("N", H_Sb_table$Var1), ]
  H_Sb_table[,3] <- H_Sb_table[,2]/sum(H_Sb_table[,2])
  colnames(H_Sb_table) <- c("base1234","H_number","H_proportion")
  
  # 4. softclip类型组合统计
  H11 <- subset(H, H$R1S5type != 0 & H$R2S5type != 0)
  numberH11 <- length(H11[,1])
  H10 <- subset(H, H$R1S5type != 0 & H$R2S5type == 0)
  numberH10 <- length(H10[,1])
  H01 <- subset(H, H$R1S5type == 0 & H$R2S5type != 0)
  numberH01 <- length(H01[,1])
  H_type3_table <- rbind(numberH11, numberH10, numberH01)
  H_type3_table <- as.data.frame(H_type3_table)
  H_type3_table[,2] <- H_type3_table[,1]/sum(H_type3_table[,1])
  H_type3_table <- cbind(row_names = rownames(H_type3_table), H_type3_table)
  colnames(H_type3_table) <- c("Htype3","Hnumber","Hproportion")
  
  # 输出文件名前缀
  file_prefix <- sub("_Hwenjian.txt", "", ff1[i])
  
  # 保存结果到对应目录
  write.table(H_Snumber_table, 
              file.path(output_base, "H_Snumber_table", paste0(file_prefix, "_H_Snumber_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(H_Sbase_table,
              file.path(output_base, "H_Sbase_table", paste0(file_prefix, "_H_Sbase_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(H_Sb_table,
              file.path(output_base, "H_Sb_table", paste0(file_prefix, "_H_Sb_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(H_type3_table,
              file.path(output_base, "H_type3_table", paste0(file_prefix, "_H_type3_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}






# 定义输入和输出路径
input_dir <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Lwenjian/"
output_base <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/softclipseq/"

# 创建输出目录（在循环外创建，避免重复创建）
for(dir_name in c("L_Snumber_table", "L_Sbase_table", "L_Sb_table", "L_type3_table")) {
  dir.create(file.path(output_base, dir_name), recursive = TRUE, showWarnings = FALSE)
}

# 获取文件列表
ff1 <- list.files(input_dir, pattern = "_Lwenjian.txt$")

# 处理每个文件
for(i in 1:length(ff1)) {
  # 读取数据
  L <- read.table(paste(input_dir, ff1[i], sep = ""), sep = "\t", head = TRUE)
  num_L <- nrow(L)
  
  # 1. softclip长度分布统计
  L_Snumber_tab <- as.data.frame(table(cut(na.omit(c(L$R1S5number, L$R2S5number)), 
                                           breaks = seq(0, 30, by = 1), 
                                           right = FALSE, include.lowest = TRUE)))
  L_Snumber_table <- L_Snumber_tab[-1,]
  L_Snumber_table[,3] <- L_Snumber_table[,2] / sum(L_Snumber_table[,2])
  colnames(L_Snumber_table) <- c("Type","Number","Proportion")
  
  # 2. softclip序列碱基组成统计
  L_Sbase <- c(L$R1S5base, L$R2S5base)
  L_Sbase_seq <- paste(L_Sbase, collapse = "")
  L_Sbase_table <- table(strsplit(L_Sbase_seq, "")[[1]])
  L_Sbase_table <- as.data.frame(L_Sbase_table)
  L_Sbase_table <- L_Sbase_table[!grepl("N", L_Sbase_table$Var1), ]
  L_Sbase_table[,3] <- L_Sbase_table[,2]/sum(L_Sbase_table[,2])
  colnames(L_Sbase_table) <- c("base1234","L_number","L_proportion")
  
  # 3. softclip第一个碱基组成统计
  L_R11 <- subset(L, L$R1S5type != 0)
  L_R21 <- subset(L, L$R2S5type != 0)
  L_Sb <- c(L_R11$R1S5b1, L_R21$R2S5b1)
  L_Sb_seq <- paste(L_Sb, collapse = "")
  L_Sb_table <- table(strsplit(L_Sb_seq, "")[[1]])
  L_Sb_table <- as.data.frame(L_Sb_table)
  L_Sb_table <- L_Sb_table[!grepl("N", L_Sb_table$Var1), ]
  L_Sb_table[,3] <- L_Sb_table[,2]/sum(L_Sb_table[,2])
  colnames(L_Sb_table) <- c("base1234","L_number","L_proportion")
  
  # 4. softclip类型组合统计
  L11 <- subset(L, L$R1S5type != 0 & L$R2S5type != 0)
  numberL11 <- length(L11[,1])
  L10 <- subset(L, L$R1S5type != 0 & L$R2S5type == 0)
  numberL10 <- length(L10[,1])
  L01 <- subset(L, L$R1S5type == 0 & L$R2S5type != 0)
  numberL01 <- length(L01[,1])
  L_type3_table <- rbind(numberL11, numberL10, numberL01)
  L_type3_table <- as.data.frame(L_type3_table)
  L_type3_table[,2] <- L_type3_table[,1]/sum(L_type3_table[,1])
  L_type3_table <- cbind(row_names = rownames(L_type3_table), L_type3_table)
  colnames(L_type3_table) <- c("Ltype3","Lnumber","Lproportion")
  
  # 输出文件名前缀
  file_prefix <- sub("_Lwenjian.txt", "", ff1[i])
  
  # 保存结果到对应目录
  write.table(L_Snumber_table, 
              file.path(output_base, "L_Snumber_table", paste0(file_prefix, "_L_Snumber_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(L_Sbase_table,
              file.path(output_base, "L_Sbase_table", paste0(file_prefix, "_L_Sbase_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(L_Sb_table,
              file.path(output_base, "L_Sb_table", paste0(file_prefix, "_L_Sb_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(L_type3_table,
              file.path(output_base, "L_type3_table", paste0(file_prefix, "_L_type3_table.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}






