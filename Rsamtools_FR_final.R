# ============================================================================
# Script: Rsamtools_FR_final.R
# Description:
#   This script processes paired-end BAM files (aligned to mitochondrial genome)
#   to extract soft-clip information, compute quality metrics, and generate
#   two tab-delimited output files per input BAM:
#     *_Hwenjian.txt: for reads with flags 99 (R1) and 147 (R2)
#     *_Lwenjian.txt: for reads with flags 163 (R2) and 83 (R1)
#   The outputs contain soft-clip positions, sequences, base qualities, GC
#   content, and fragment start/end positions.
# Input:
#   "_finally.mt.bam" and "_finally.mt.bam.bai"
# Output:
#   "_Hwenjian.txt" and "_Lwenjian.txt"
# Usage: Rscript Rsamtools_FR_final.R
# ============================================================================




setwd("D:/run/mtbam_mtbambai/")
library(Rsamtools)
library(stringr)


# List BAM files and their index files
bam<-list.files("D:/run/mtbam_mtbambai/",pattern = "\\_finally.mt.bam$")
bam
index<-list.files("D:/run/mtbam_mtbambai/",pattern = "\\_finally.mt.bam.bai$")
index


# Process each BAM file
for(i in 1:length(bam)){

bamFile <- bam[i]
indexFile <- index[i]
bamData <- scanBam(bamFile)
A<-as.data.frame(bamData[[1]])
head(A)
dim(A)


# Absolute insert size
A$isize <- abs(A$isize)


# Extract 5' soft-clip length from CIGAR string
A[,14]<-str_extract(A$cigar,"^[0-9]+S")
A[,15]<-str_extract(A[,14],"^[0-9]+")
A[,15]<-as.numeric(A[,15])
A$V14[is.na(A$V14)] <- 0
A$V15[is.na(A$V15)] <- 0
# Extract soft-clip bases at 5' end
A[,16]<-substr(A$seq, 1, A[,15])               # soft-clip region
A[,17]<-substr(A$seq, A[,15]+1, A[,15]+1)      # first base after soft-clip
A[,18]<-substr(A$seq, A[,15]+1, A[,15]+2)      # first two bases after soft-clip
A[,19]<-substr(A$seq, A[,15]+1, A[,15]+3)      # first three bases after soft-clip
A[,20]<-substr(A$seq, A[,15]+1, A[,15]+4)      # first four bases after soft-clip


# Extract 3' soft-clip length from CIGAR string
A[,21]<-str_extract(A$cigar,"[0-9]+S$")
A[,22]<-str_extract(A[,21],"^[0-9]+")
A[,22]<-as.numeric(A[,22])
A$V21[is.na(A$V21)] <- 0
A$V22[is.na(A$V22)] <- 0
# Extract soft-clip bases at 3' end
A[,23]<-substr(A$seq, nchar(A$seq)-A[,22]+1, nchar(A$seq))             # soft-clip region
A[,24]<-substr(A$seq, nchar(A$seq)-A[,22], nchar(A$seq)-A[,22])        # first base before soft-clip
A[,25]<-substr(A$seq, nchar(A$seq)-A[,22]-1, nchar(A$seq)-A[,22])      # first two bases before soft-clip
A[,26]<-substr(A$seq, nchar(A$seq)-A[,22]-2, nchar(A$seq)-A[,22])      # first three bases before soft-clip
A[,27]<-substr(A$seq, nchar(A$seq)-A[,22]-3, nchar(A$seq)-A[,22])      # first four bases before soft-clip




# Function to extract and sum all numbers from a string (e.g., CIGAR)
extract_and_sum_numbers <- function(str) {
  numbers_as_strings <- unlist(strsplit(str, "[^0-9]+"))
  sum_of_numbers <- sum(as.numeric(numbers_as_strings))
  return(sum_of_numbers)
}
sums <- sapply(A$cigar, extract_and_sum_numbers)
sums_copy<-as.data.frame(sums)


# Total aligned length, compute end position
A[,28]<-sums_copy$sums
A[,29]<-A[,5]+A[,28]-1-A[,15]-A[,22]


# Split quality strings into segments: 5' soft-clip, matched, 3' soft-clip
A[,30]<-substr(A$qual, 1, A[,15])                                  # 5' soft-clip
A[,31]<-substr(A$qual, A[,15]+1, nchar(A$qual)-A[,22])             # matched
A[,32]<-substr(A$qual, nchar(A$qual)-A[,22]+1, nchar(A$qual))      # 3' soft-clip




# Convert ASCII quality to Phred scores
qual_to_phred <- function(qual_str) {
  if (is.na(qual_str)) return(NA_integer_)
  as.integer(charToRaw(qual_str)) - 33
}
A$V30 <- lapply(A$V30, qual_to_phred)
A$V31 <- lapply(A$V31, qual_to_phred)
A$V32 <- lapply(A$V32, qual_to_phred)




# Compute mean Phred scores for each segment
A$V33 <- sapply(A$V30, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})
A$V34 <- sapply(A$V31, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})
A$V35 <- sapply(A$V32, function(x) {
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
})
A$V33[is.na(A$V33)] <- 0
A$V34[is.na(A$V34)] <- 0
A$V35[is.na(A$V35)] <- 0

# Get the matched portion of the sequence (excluding soft-clips)
A[,36]<-substr(A$seq, A[,15]+1, nchar(A$seq)-A[,22])

# Compute GC percentage of the matched portion
A[,37] <- round((str_count(A$V36, "G") + str_count(A$V36, "C")) / nchar(A$V36) * 100, 2)




# ------------------- Process H-strand pairs (flags 99 and 147) -------------------
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

# Combine R1 (flag 99) and R2 (flag 147) by qname
H_new<-cbind(H_99copy[,c(1:11,14:17,21:24,28:35,37)],H_147copy[,c(1:11,14:17,21:24,28:35,37)])
names(H_new) <- c("R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S5type","R1S5number","R1S5base","R1S5b1","R1S3type","R1S3number","R1S3base","R1S3b1","R1length","R1endpos","R1S5qual","R1Mqual","R1S3qual","R1S5qmean","R1Mqmean","R1S3qmean","R1GC","R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S3type","R2S3number","R2S3base","R2S3b1","R2S5type","R2S5number","R2S5base","R2S5b1","R2length","R2endpos","R2S3qual","R2Mqual","R2S5qual","R2S3qmean","R2Mqmean","R2S5qmean","R2GC")

# Compute fragment start and end positions
H_new[,57]<-apply(H_new[,c(5,21,33,49)],1,min)
H_new[,58]<-apply(H_new[,c(5,21,33,49)],1,max)
names(H_new) <- c("R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S5type","R1S5number","R1S5base","R1S5b1","R1S3type","R1S3number","R1S3base","R1S3b1","R1length","R1endpos","R1S5qual","R1Mqual","R1S3qual","R1S5qmean","R1Mqmean","R1S3qmean","R1GC","R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S3type","R2S3number","R2S3base","R2S3b1","R2S5type","R2S5number","R2S5base","R2S5b1","R2length","R2endpos","R2S3qual","R2Mqual","R2S5qual","R2S3qmean","R2Mqmean","R2S5qmean","R2GC","startpos","endpos")

H_new<-H_new[order(H_new$startpos),]
H_final<-H_new[,c(1,57,58,11:19,25:28,39:47,53:56)]

# Create output directory for H-strand files and write file
dir.create("D:/run/Hwenjian/")
write.table(x = H_final, file = paste0("D:/run/Hwenjian/", sub("_finally\\.mt\\.bam$", "_Hwenjian.txt", bam[i])), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




# ------------------- Process L-strand pairs (flags 163 and 83) -------------------
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

# Combine R2 (flag 163) and R1 (flag 83) by qname
L_new<-cbind(L_163copy[,c(1:11,14:17,21:24,28:35,37)],L_83copy[,c(1:11,14:17,21:24,28:35,37)])
names(L_new) <- c("R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S5type","R2S5number","R2S5base","R2S5b1","R2S3type","R2S3number","R2S3base","R2S3b1","R2length","R2endpos","R2S5qual","R2Mqual","R2S3qual","R2S5qmean","R2Mqmean","R2S3qmean","R2GC","R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S3type","R1S3number","R1S3base","R1S3b1","R1S5type","R1S5number","R1S5base","R1S5b1","R1length","R1endpos","R1S3qual","R1Mqual","R1S5qual","R1S3qmean","R1Mqmean","R1S5qmean","R1GC")

# Compute fragment start and end positions
L_new[,57]<-apply(L_new[,c(5,21,33,49)],1,min)
L_new[,58]<-apply(L_new[,c(5,21,33,49)],1,max)
names(L_new) <- c("R2qname","R2flag","R2rname","R2strand","R2pos","R2qwidth","R2mapq","R2cigar","R2mrnm","R2mpos","R2isize","R2S5type","R2S5number","R2S5base","R2S5b1","R2S3type","R2S3number","R2S3base","R2S3b1","R2length","R2endpos","R2S5qual","R2Mqual","R2S3qual","R2S5qmean","R2Mqmean","R2S3qmean","R2GC","R1qname","R1flag","R1rname","R1strand","R1pos","R1qwidth","R1mapq","R1cigar","R1mrnm","R1mpos","R1isize","R1S3type","R1S3number","R1S3base","R1S3b1","R1S5type","R1S5number","R1S5base","R1S5b1","R1length","R1endpos","R1S3qual","R1Mqual","R1S5qual","R1S3qmean","R1Mqmean","R1S5qmean","R1GC","startpos","endpos")

L_new<-L_new[order(L_new$startpos),]
L_final<-L_new[,c(29,57,58,11:19,25:28,39:47,53:56)]

# Create output directory for L-strand files and write file
dir.create("D:/run/Lwenjian/")
write.table(x = L_final, file = paste0("D:/run/Lwenjian/", sub("_finally\\.mt\\.bam$", "_Lwenjian.txt", bam[i])), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}



