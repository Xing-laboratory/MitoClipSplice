# ====================== 1. 统一设置路径变量（核心修改） ======================
# 根路径：所有路径的基础，修改时只需改这一处
root_path <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/"
# 输入文件子目录
input_subdir <- "Hwenjian/"
# 完整输入路径（拼接根路径+输入子目录）
input_full_path <- paste0(root_path, input_subdir)
# 输出父子目录（新增的softclipfrag层）
output_parent_subdir <- "softclipfrag/"
# 完整输出父路径
output_root_path <- paste0(root_path, output_parent_subdir)

# 定义所有输出子目录（便于统一管理和创建）
output_subdirs <- c(
  "H1_Mqmean_table",
  "H0_Mqmean_table",
  "H1_GC_table",
  "H0_GC_table",
  "H1_insert_table",
  "H0_insert_table",
  "H_result"
)

# ====================== 2. 提前创建所有输出目录（避免循环内重复创建） ======================
# 循环创建所有输出子目录，showWarnings=FALSE避免目录已存在时的警告
for (subdir in output_subdirs) {
  # 拼接完整输出子目录路径
  dir_path <- paste0(output_root_path, subdir, "/")
  # recursive=TRUE确保嵌套目录（如softclipfrag）也能被创建
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ====================== 3. 核心分析逻辑（使用路径变量） ======================
# 获取输入目录下的目标文件
ff1 <- list.files(input_full_path, pattern = "\\_Hwenjian.txt$")
cat("找到待处理的文件数量：", length(ff1), "\n")
print(ff1)

# 循环处理每个文件
for(i in 1:length(ff1)){
  # 读取文件（使用路径变量拼接，替代硬编码）
  H <- read.table(paste0(input_full_path, ff1[i]), sep = "\t", header = TRUE)
  # 更推荐用nrow()获取数据框行数（比length(H[,1])更直观、健壮）
  num_H <- nrow(H)
  
  # 按soft-clipped分类筛选
  H1 <- subset(H, H$R1S5type != 0 | H$R2S5type != 0)
  number_H1 <- nrow(H1)
  number_H1_pro <- number_H1 / num_H
  
  H0 <- subset(H, H$R1S5type == 0 & H$R2S5type == 0)
  number_H0 <- nrow(H0)
  number_H0_pro <- number_H0 / num_H
  
  # -------------------- Mqmean 分析 --------------------
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
  
  # 计算Q30比例（30-50区间，对应breaks的16-25行）
  H1_Mqmean_pro <- sum(H1_Mqmean_table[31:50,3]) / sum(H1_Mqmean_table[,3])
  H0_Mqmean_pro <- sum(H0_Mqmean_table[31:50,3]) / sum(H0_Mqmean_table[,3])
  
  # 写入Mqmean结果（使用路径变量拼接）
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
  
  # -------------------- GC 分析 --------------------
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
  
  # 计算GC中位数（增加na.rm=TRUE避免NA值导致结果出错）
  H1_GC_median <- median(c(H1$R1GC, H1$R2GC), na.rm = TRUE)
  H0_GC_median <- median(c(H0$R1GC, H0$R2GC), na.rm = TRUE)
  
  # 写入GC结果
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
  
  # -------------------- Insert size 分析 --------------------
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
  
  # 计算insert size中位数（增加na.rm=TRUE）
  H1_insert_median <- median(H1$R1isize, na.rm = TRUE)
  H0_insert_median <- median(H0$R1isize, na.rm = TRUE)
  
  # 写入insert size结果
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
  
  # -------------------- 汇总结果 --------------------
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
  
  # 写入汇总结果
  write.table(
    x = H_result, 
    file = paste0(output_root_path, "H_result/", sub("_Hwenjian.txt", "_H_result.txt", ff1[i])), 
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
  )
  
  # 进度提示：打印当前处理的文件，便于监控运行状态
  cat(paste("已完成：", i, "/", length(ff1), " 文件：", ff1[i], "\n"))
}

cat("所有文件处理完成！结果已输出到：", output_root_path, "\n")






# ====================== 1. 统一设置路径变量（核心修改） ======================
# 根路径：所有路径的基础，修改时只需改这一处
root_path <- "D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/"
# 输入文件子目录
input_subdir <- "Lwenjian/"
# 完整输入路径（拼接根路径+输入子目录）
input_full_path <- paste0(root_path, input_subdir)
# 输出父子目录（新增的softclipfrag层）
output_parent_subdir <- "softclipfrag/"
# 完整输出父路径
output_root_path <- paste0(root_path, output_parent_subdir)

# 定义所有输出子目录（便于统一管理和创建）
output_subdirs <- c(
  "L1_Mqmean_table",
  "L0_Mqmean_table",
  "L1_GC_table",
  "L0_GC_table",
  "L1_insert_table",
  "L0_insert_table",
  "L_result"
)

# ====================== 2. 提前创建所有输出目录（避免循环内重复创建） ======================
# 循环创建所有输出子目录，showWarnings=FALSE避免目录已存在时的警告
for (subdir in output_subdirs) {
  # 拼接完整输出子目录路径
  dir_path <- paste0(output_root_path, subdir, "/")
  # recursive=TRUE确保嵌套目录（如softclipfrag）也能被创建
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ====================== 3. 核心分析逻辑（使用路径变量） ======================
# 获取输入目录下的目标文件
ff1 <- list.files(input_full_path, pattern = "\\_Lwenjian.txt$")
cat("找到待处理的文件数量：", length(ff1), "\n")
print(ff1)

# 循环处理每个文件
for(i in 1:length(ff1)){
  # 读取文件（使用路径变量拼接，替代硬编码）
  L <- read.table(paste0(input_full_path, ff1[i]), sep = "\t", header = TRUE)
  # 更推荐用nrow()获取数据框行数（比length(L[,1])更直观、健壮）
  num_L <- nrow(L)
  
  # 按soft-clipped分类筛选
  L1 <- subset(L, L$R1S5type != 0 | L$R2S5type != 0)
  number_L1 <- nrow(L1)
  number_L1_pro <- number_L1 / num_L
  
  L0 <- subset(L, L$R1S5type == 0 & L$R2S5type == 0)
  number_L0 <- nrow(L0)
  number_L0_pro <- number_L0 / num_L
  
  # -------------------- Mqmean 分析 --------------------
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
  
  # 计算Q30比例（30-50区间，对应breaks的16-25行）
  L1_Mqmean_pro <- sum(L1_Mqmean_table[31:50,3]) / sum(L1_Mqmean_table[,3])
  L0_Mqmean_pro <- sum(L0_Mqmean_table[31:50,3]) / sum(L0_Mqmean_table[,3])
  
  # 写入Mqmean结果（使用路径变量拼接）
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
  
  # -------------------- GC 分析 --------------------
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
  
  # 计算GC中位数（增加na.rm=TRUE避免NA值导致结果出错）
  L1_GC_median <- median(c(L1$R1GC, L1$R2GC), na.rm = TRUE)
  L0_GC_median <- median(c(L0$R1GC, L0$R2GC), na.rm = TRUE)
  
  # 写入GC结果
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
  
  # -------------------- Insert size 分析 --------------------
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
  
  # 计算insert size中位数（增加na.rm=TRUE）
  L1_insert_median <- median(L1$R1isize, na.rm = TRUE)
  L0_insert_median <- median(L0$R1isize, na.rm = TRUE)
  
  # 写入insert size结果
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
  
  # -------------------- 汇总结果 --------------------
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
  
  # 写入汇总结果
  write.table(
    x = L_result, 
    file = paste0(output_root_path, "L_result/", sub("_Lwenjian.txt", "_L_result.txt", ff1[i])), 
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
  )
  
  # 进度提示：打印当前处理的文件，便于监控运行状态
  cat(paste("已完成：", i, "/", length(ff1), " 文件：", ff1[i], "\n"))
}

cat("所有文件处理完成！结果已输出到：", output_root_path, "\n")






