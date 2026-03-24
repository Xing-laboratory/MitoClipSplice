setwd("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/mtbam_mtbambai/")
library(Rsamtools)
library(stringr)

bam<-list.files("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/mtbam_mtbambai/",pattern = "\\_finally.mt.bam$")
bam

index<-list.files("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/mtbam_mtbambai/",pattern = "\\_finally.mt.bam.bai$")
index


for(i in 1:length(bam)){

bamFile <- bam[i]
indexFile <- index[i]
bamData <- scanBam(bamFile)
A<-as.data.frame(bamData[[1]])
head(A)
dim(A)

A$isize <- abs(A$isize)

A[,14]<-str_extract(A$cigar,"^[0-9]+S")
A[,15]<-str_extract(A[,14],"^[0-9]+")
A[,15]<-as.numeric(A[,15])
A$V14[is.na(A$V14)] <- 0
A$V15[is.na(A$V15)] <- 0
A[,16]<-substr(A$seq, 1, A[,15])
A[,17]<-substr(A$seq, A[,15]+1, A[,15]+1)
A[,18]<-substr(A$seq, A[,15]+1, A[,15]+2)
A[,19]<-substr(A$seq, A[,15]+1, A[,15]+3)
A[,20]<-substr(A$seq, A[,15]+1, A[,15]+4)


A[,21]<-str_extract(A$cigar,"[0-9]+S$")
A[,22]<-str_extract(A[,21],"^[0-9]+")
A[,22]<-as.numeric(A[,22])
A$V21[is.na(A$V21)] <- 0
A$V22[is.na(A$V22)] <- 0
A[,23]<-substr(A$seq, nchar(A$seq)-A[,22]+1, nchar(A$seq))
A[,24]<-substr(A$seq, nchar(A$seq)-A[,22], nchar(A$seq)-A[,22])
A[,25]<-substr(A$seq, nchar(A$seq)-A[,22]-1, nchar(A$seq)-A[,22])
A[,26]<-substr(A$seq, nchar(A$seq)-A[,22]-2, nchar(A$seq)-A[,22])
A[,27]<-substr(A$seq, nchar(A$seq)-A[,22]-3, nchar(A$seq)-A[,22])




extract_and_sum_numbers <- function(str) {  
  # 使用正则表达式提取所有数字作为字符串  
  numbers_as_strings <- unlist(strsplit(str, "[^0-9]+"))  
  # 将字符串转换为数值并求和  
  sum_of_numbers <- sum(as.numeric(numbers_as_strings))  
  # 返回求和结果  
  return(sum_of_numbers)  
}  

# 使用sapply对每个字符串应用上述函数  
sums <- sapply(A$cigar, extract_and_sum_numbers) 
sums_copy<-as.data.frame(sums)

A[,28]<-sums_copy$sums
A[,29]<-A[,5]+A[,28]-1-A[,15]-A[,22]


A[,30]<-substr(A$qual, 1, A[,15])
A[,31]<-substr(A$qual, A[,15]+1, nchar(A$qual)-A[,22])
A[,32]<-substr(A$qual, nchar(A$qual)-A[,22]+1, nchar(A$qual))




# 转换函数：单个qual字符串 → Phred值向量（处理NA）
qual_to_phred <- function(qual_str) {
  if (is.na(qual_str)) return(NA_integer_)
  as.integer(charToRaw(qual_str)) - 33  # Phred+33编码（主流）
}

# 批量转换所有read的qual字符串
A$V30 <- lapply(A$V30, qual_to_phred)
A$V31 <- lapply(A$V31, qual_to_phred)
A$V32 <- lapply(A$V32, qual_to_phred)




# 计算每个read的平均Phred值（排除NA）
A$V33 <- sapply(A$V30, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})


# 计算每个read的平均Phred值（排除NA）
A$V34 <- sapply(A$V31, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})


# 计算每个read的平均Phred值（排除NA）
A$V35 <- sapply(A$V32, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})


A$V33[is.na(A$V33)] <- 0
A$V34[is.na(A$V34)] <- 0
A$V35[is.na(A$V35)] <- 0

A[,36]<-substr(A$seq, A[,15]+1, nchar(A$seq)-A[,22])

# 直接新增gc_ratio列，无中间列
A[,37] <- round((str_count(A$V36, "G") + str_count(A$V36, "C")) / nchar(A$V36) * 100, 2)




H <- subset(A,A$flag == 99 | A$flag == 147)
dim(H)
H_99 <- subset(A,A$flag == 99)
dim(H_99)
H_147 <- subset(A,A$flag == 147)
dim(H_147)

H_99copy<-H_99[order(H_99$qname),]
dim(H_99copy)
H_147copy<-H_147[order(H_147$qname),]
dim(H_147copy)

H_new<-cbind(H_99copy[,c(1:11,14:17,21:24,28:35,37)],H_147copy[,c(1:11,14:17,21:24,28:35,37)])
names(H_new) <- c("R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S5type","R1S5number","R1S5base","R1S5b1","R1S3type","R1S3number","R1S3base","R1S3b1","R1length","R1endpos","R1S5qual","R1Mqual","R1S3qual","R1S5qmean","R1Mqmean","R1S3qmean","R1GC","R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S3type","R2S3number","R2S3base","R2S3b1","R2S5type","R2S5number","R2S5base","R2S5b1","R2length","R2endpos","R2S3qual","R2Mqual","R2S5qual","R2S3qmean","R2Mqmean","R2S5qmean","R2GC")

H_new[,57]<-apply(H_new[,c(5,21,33,49)],1,min)
H_new[,58]<-apply(H_new[,c(5,21,33,49)],1,max)
names(H_new) <- c("R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S5type","R1S5number","R1S5base","R1S5b1","R1S3type","R1S3number","R1S3base","R1S3b1","R1length","R1endpos","R1S5qual","R1Mqual","R1S3qual","R1S5qmean","R1Mqmean","R1S3qmean","R1GC","R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S3type","R2S3number","R2S3base","R2S3b1","R2S5type","R2S5number","R2S5base","R2S5b1","R2length","R2endpos","R2S3qual","R2Mqual","R2S5qual","R2S3qmean","R2Mqmean","R2S5qmean","R2GC","startpos","endpos")

H_new<-H_new[order(H_new$startpos),]
H_final<-H_new[,c(1,57,58,11:19,25:28,39:47,53:56)]
dir.create("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Hwenjian/")
# 核心：将目录、替换后的文件名拼接为完整的文件路径，作为write.table的第二个参数
write.table(x = H_final, file = paste0("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Hwenjian/", sub("_finally\\.mt\\.bam$", "_Hwenjian.txt", bam[i])), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




L <- subset(A,A$flag == 163 | A$flag == 83)
dim(L)
L_163 <- subset(A,A$flag == 163)
dim(L_163)
L_83 <- subset(A,A$flag == 83)
dim(L_83)

L_163copy<-L_163[order(L_163$qname),]
dim(L_163copy)
L_83copy<-L_83[order(L_83$qname),]
dim(L_83copy)

L_new<-cbind(L_163copy[,c(1:11,14:17,21:24,28:35,37)],L_83copy[,c(1:11,14:17,21:24,28:35,37)])
names(L_new) <- c("R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S5type","R2S5number","R2S5base","R2S5b1","R2S3type","R2S3number","R2S3base","R2S3b1","R2length","R2endpos","R2S5qual","R2Mqual","R2S3qual","R2S5qmean","R2Mqmean","R2S3qmean","R2GC","R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S3type","R1S3number","R1S3base","R1S3b1","R1S5type","R1S5number","R1S5base","R1S5b1","R1length","R1endpos","R1S3qual","R1Mqual","R1S5qual","R1S3qmean","R1Mqmean","R1S5qmean","R1GC")

L_new[,57]<-apply(L_new[,c(5,21,33,49)],1,min)
L_new[,58]<-apply(L_new[,c(5,21,33,49)],1,max)
names(L_new) <- c("R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S5type","R2S5number","R2S5base","R2S5b1","R2S3type","R2S3number","R2S3base","R2S3b1","R2length","R2endpos","R2S5qual","R2Mqual","R2S3qual","R2S5qmean","R2Mqmean","R2S3qmean","R2GC","R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S3type","R1S3number","R1S3base","R1S3b1","R1S5type","R1S5number","R1S5base","R1S5b1","R1length","R1endpos","R1S3qual","R1Mqual","R1S5qual","R1S3qmean","R1Mqmean","R1S5qmean","R1GC","startpos","endpos")

L_new<-L_new[order(L_new$startpos),]
L_final<-L_new[,c(29,57,58,11:19,25:28,39:47,53:56)]
dir.create("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Lwenjian/")
write.table(x = L_final, file = paste0("D:/run/softclip_newcode_260130_网站/FR_WTS_private_random_barcode_Homo/Lwenjian/", sub("_finally\\.mt\\.bam$", "_Lwenjian.txt", bam[i])), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}




