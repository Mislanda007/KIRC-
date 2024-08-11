data<-read.csv("D:/生物医学工程/test1 (1).csv",row.names=1)
View(data)
library(readxl)
clin <- read_excel("D:/生物医学工程/TCGR-KIRC_clinical.xlsx")
View(clin)
#coad.clin <- clin[clin$OS.time>=30,]
coad.clin <- clin[!duplicated(clin$sample),]#去重
View(coad.clin)
protein <- read_xlsx( "D:/生物医学工程/immune_GeneList.xlsx")
View(protein)
dim(protein)
gene <- intersect(protein$Symbol,rownames(data))
sample<-coad.clin$sample
sample=gsub("-",".",sample)
sample <- intersect(sample,colnames(data))
enengyTurExp <- data[gene,sample]
dim(enengyTurExp)
View(enengyTurExp)
install.packages("writexl")

# 加载 writexl 包
library(writexl)

# 将数据保存为 Excel 文件
write_xlsx(data, path = "D:/生物医学工程/enengyTurExp_up2.xlsx")
library(tidyr)
genetable <- lapply(as.data.frame(t(enengyTurExp)),function(r){sum(r==0)}) %>% as.data.frame()
enengyTurExp <- enengyTurExp[colnames(genetable)[genetable[1,] < length(colnames(enengyTurExp))/2],]
sum(is.na(enengyTurExp))
enengyTurExp <- na.omit(enengyTurExp)
sum(is.na(enengyTurExp))
dim(enengyTurExp)
write_xlsx(enengyTurExp, path = "D:/生物医学工程/enengyTurExp_up2.xlsx")
library(NMF)
coad.log2fpkm.enengy <- enengyTurExp##
ranks <- 2:10
estim.coad <- nmf(coad.log2fpkm.enengy,ranks, nrun=50)
duplicated(colnames(coad.log2fpkm.enengy))
plot(estim.coad)
seed = 2020820
nmf.rank4 <- nmf(coad.log2fpkm.enengy, 
                 rank = 3, 
                 nrun=50,
                 seed = seed, 
                 method = "brunet")
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
index <- extractFeatures(nmf.rank4,"max") 
sig.order <- unlist(index)
NMF.Exp.rank4 <- coad.log2fpkm.enengy[sig.order,]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4) #sig.order有时候会有缺失值
group <- predict(nmf.rank4) # 提出亚型
table(group)
#data<-data.frame(group)
write.table(group, file  = "D:/生物医学工程/groupup2.txt",sep=" ")
consensusmap(nmf.rank4,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank4)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
et<-rownames(enengyTurExp)
write.table(et, file  = "D:/生物医学工程/id3.txt",sep=" ")
