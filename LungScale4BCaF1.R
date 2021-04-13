#

#     LungScale4BCaF1.R
#     /DS/RhDec20/LungScale4/
#      December 9, 2020
#
#   Run Dec 31 for pdfs for NS
#  modified to use all AF data in scoring
#
#  January 5, 2021 to output Lung PCA & scores
#
#  March 25, 2021 to output AF tissue scores on Lung Scale
#
#  PCA approach using results from 
#  limma fit (BCaF1) on batch-corrected data
#
library(ggplot2)
library(PCAtools)

#   set location of the directory containing the datasets
in1DataPath="/DS/RhDec20/Filtering/"
in2DataPath="/DS/RhDec20/RefinedDatasets/RemoveExpSamples/"
in3DataPath="/DS/RhDec20/LungScale4/"
in4DataPath="/DS/RhDec20/BatchCorr/AF/"
outDataPath="/DS/RhDec20/LungScale4/Results/"
setwd(in1DataPath)

#  = = = = BRING IN log2TPM & META DATA = = = = = 

f=paste(in1DataPath,"Lung_log2TPM_BCaF1.txt",sep="")
log2TPM=read.table(file=f,sep="\t",header = T,
               check.names = F,stringsAsFactors = F)

f=paste0(in2DataPath,"MetaForTPM171819Usable.txt",sep="")
Meta.all=read.table(file=f,header = T,sep="\t",
                stringsAsFactors = F)
rownames(Meta.all)=Meta.all$RNASeqID

rownames(log2TPM)=log2TPM$UID
log2TPM$UID=NULL


# = = = SPECIFY THE TISSUE & SELECT DATA = = = = =
Tissues=c("Lung","Brain","AF")
Tis=Tissues[1]

mT=Meta.all[Meta.all$Tissue==Tis,]

# -- Specify the expression data matrix
expT=log2TPM[,mT$RNASeqID]

# -- Remove 0 rows (0 rows mess up z-scoring)
expT=expT[rowSums(expT) > 0,]


# ==== BRING IN THE GENE LIST ======
Genelist="Lung4aBCaF1"  #  adjp <= 0.05 & abs(lfc) >=log2(2)
f=paste0(in3DataPath,"LungGenelist4a.rds")
Gene.list=readRDS(f)
Gene.list[["Up"]]=Gene.list[["Up"]][order(Gene.list[["Up"]])]
Gene.list[["Down"]]=Gene.list[["Down"]][order(Gene.list[["Down"]])]

Gene.list[["Up"]]
#   Up list (a) includes SFTPA1, SFTPC, CTSH
Gene.list[["Down"]]
#  Downlist (a) includes 7 collagens & some cell cycle

#  Possible some features not in expT due to filtering
#  Remove those from lists
GL=list()
GL[["Up"]]=Gene.list[["Up"]][(Gene.list[["Up"]] %in% rownames(expT))]
GL[["Down"]]=Gene.list[["Down"]][(Gene.list[["Down"]] %in% rownames(expT))]
#  all present
# 
expT.up=expT[GL[["Up"]],]
expT.up=expT.up[order(rownames(expT.up)),]
expT.down=expT[GL[["Down"]],]
expT.down=expT.down[order(rownames(expT.down)),]


# ====== PCA Approach =================

#  Do PCA on combined set
expT.both=rbind.data.frame(expT.up,expT.down,
               stringsAsFactors = F)

# -- Apply only to TC & PTC data
mT.TCPTC=mT[mT$TrtShort %in% c("PTC","TC"), ]
rownames(mT.TCPTC)=mT.TCPTC$RNASeqID
expT.both.TCPTC=expT.both[,mT.TCPTC$RNASeqID]

pca1=pca(expT.both.TCPTC,
         metadata = mT.TCPTC,
         removeVar = 0,scale=T,center=T)
tl=paste0(Tis," (BCaF1) data, TC & PTC samples, PCA-based scoring")
st=paste0("Genelist used:",Genelist)
screeplot(pca1,title = tl,subtitle = st)
biplot(pca1,colby = 'TrtShort',shape = 'Sex',
       legendPosition="top",title=tl,
       subtitle = st)

pca.obj=pca1
pc1.loadings=pca.obj$loadings$PC1

#  Score all the samples on pc1
exp.z=t(scale(t(expT),center = T,scale=T))
summary(rowMeans(exp.z))
summary(apply(exp.z,1,sd))

#  From PCAtools vignette
p.prcomp <- list(sdev = pca.obj$sdev, rotation = data.matrix(pca.obj$loadings), 
                 x = data.matrix(pca.obj$rotated),
                    center = TRUE, scale = TRUE)
class(p.prcomp) <- 'prcomp'
#  provide entire Lung dataset
newdata <- t(exp.z) 
Preds=predict(p.prcomp, newdata = newdata)[,1:5]   

MB=merge(mT,Preds,
         by='row.names') 
MB$TrtShortF=factor(MB$TrtShort)

# ***  SAVE Lung pca.obj & MB for other work ***
f=paste0(outDataPath,"LungScalePCAstuff.RData")
TCvPTC_pca.obj=pca.obj
ScoreAllSamples=MB
save(file = f,TCvPTC_pca.obj,Gene.list,
     ScoreAllSamples)
# ****************************************

tl=paste0(Tis," (BCaF1) data, all samples, PCA-based scoring")
st=paste0("Genelist used:",Genelist)
href=median(MB$PC1[MB$TrtShort=="PTC"])
scoring.bp=ggplot(MB,aes(x=TrtShortF,y=PC1)) +
   geom_boxplot(outlier.shape = NA) +
   geom_jitter(width=0.2,height=0,aes(color=Sex)) +
   labs(title=tl,subtitle = st) +
   geom_hline(yintercept = href,color="black",
              lty="dotted")
scoring.bp

#****************************************
# = = = SCORE AF SAMPLES ON Lung Scale = = = 
#  !!! Use unfiltered AF data
# ***************************************
Tissues=c("Lung","Brain","AF")
Tis=Tissues[3]

#  Extract AF META DATA
mT=Meta.all[Meta.all$Tissue==Tis,]
rownames(mT)=mT$RNASeqID

# BRING IN THE AF DATA -- unfiltered
f=paste0(in4DataPath,"AF_log2TPM_BCa.txt")
AF_BC=read.table(file=f,sep="\t",header=T,
       check.names = F,stringsAsFactors = F)
rownames(AF_BC)=AF_BC$UID
AF_BC$UID=NULL
expT=AF_BC

#  Possible some features not in expT due to all 0s in AF
#  Remove those from lists
GL[["Up"]]=Gene.list[["Up"]][(Gene.list[["Up"]] %in% rownames(expT))]
GL[["Down"]]=Gene.list[["Down"]][(Gene.list[["Down"]] %in% rownames(expT))]
# lost 8

expT.up=expT[GL[["Up"]],]
expT.up=expT.up[order(rownames(expT.up)),]
expT.down=expT[GL[["Down"]],]
expT.down=expT.down[order(rownames(expT.down)),]

#  combined set
expT.both=rbind.data.frame(expT.up,expT.down,
                           stringsAsFactors = F)
exp.z=t(scale(t(expT.both),center = T,scale=T))
sum(is.nan(exp.z[,1]))

#  FIXUP FOR GENES NOT IN AF
Lost.genes=rownames(expT.both.TCPTC)
Lost.genes=Lost.genes[!(Lost.genes %in% rownames(exp.z))]
Lost.rows=matrix(data=rnorm(ncol(exp.z)*length(Lost.genes),0,0.001),
                   ncol=ncol(exp.z))
rownames(Lost.rows)=Lost.genes
colnames(Lost.rows)=colnames(exp.z)
exp.z=rbind.data.frame(exp.z,Lost.rows,stringsAsFactors = F)
exp.z=exp.z[order(rownames(exp.z)),]

summary(rowMeans(exp.z))
summary(apply(exp.z,1,sd))
newdata <- t(exp.z) 
Preds.AF=predict(p.prcomp, newdata = newdata)[,1:5]   

MB=merge(mT,Preds.AF,
         by='row.names') 
MB$TrtShortF=factor(MB$TrtShort)

tl=paste0(Tis," data (BCaF0), all samples, PCA-based scoring")
st=paste0("Genelist used:",Genelist)
href=median(MB$PC1[MB$TrtShort=="PTC"])
scoring.bp=ggplot(MB,aes(x=TrtShortF,y=PC1)) +
   geom_boxplot(outlier.shape = NA) +
   geom_jitter(width=0.2,height=0,aes(color=Sex)) +
   labs(title=tl,subtitle = st) +
   geom_hline(yintercept = href,color="black",
              lty="dotted")
scoring.bp
#  Not much signal carried over to AF

# ***  SAVE AF pca.obj & MB for other work ***
f=paste0(outDataPath,"AFonLungScalePCAstuff.RData")
AFonLung_TCvPTC_pca.obj=pca.obj
AFonLung_ScoreAllSamples=MB
save(file = f,AFonLung_TCvPTC_pca.obj,Gene.list,
     AFonLung_ScoreAllSamples)
# ****************************************

1
2
3
quit()
