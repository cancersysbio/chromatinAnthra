# compute doxo signatures from Heiser dataset (input for figure S2)

library(limma)
responses_raw = read.table("../data/heiser/procesedResponseData.tsv.csv",header = T,sep = "\t",stringsAsFactors = F)
array_exp = read.table("../data/heiser/DREAM7_DrugSensitivity1_GeneExpression.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T)

responses_raw$validNames = make.names(responses_raw$CellLine)

cellLines = colnames(array_exp[2:47])


qt = quantile(responses_raw$Doxo,c(1/3,2/3),na.rm=T)
group = as.factor(ifelse(responses_raw$Doxo>qt[2],"sens",ifelse(responses_raw$Doxo<qt[1],"res",NA)))
cellLines_doxo = responses_raw$validNames[!is.na(group)]


z = intersect(cellLines_doxo,colnames(array_exp))
idx = match(z,cellLines_doxo)
idy = match(z,colnames(array_exp))

design = model.matrix(~0+group[!is.na(group)][idx])
colnames(design) =c("res","sens")

targets = data.frame(Sample_Name = colnames(array_exp)[idy],Sample_Group = group[!is.na(group)][idx],stringsAsFactors = F)

aWeights <- arrayWeights(array_exp[,idy],design=design)
fit = lmFit(array_exp[,idy], design, weights=aWeights)
contrast.matrix = makeContrasts(res-sens,levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

fit2$genes = array_exp$HGNC_ID
sig2 <- topTable(fit2, adjust.method="BH", coef=1,number = Inf,resort.by="p")  # ESTE
sig2_noSort = topTable(fit2, adjust.method="BH", coef=1,number = Inf,sort.by="none",resort.by=NULL)  # ESTE


# generate null model
limmaNull <- function(copy, targets, nPermutations=1000){
  
  #Estimates relative quality weights for each array in a multi-array experiment.
  groups <- as.factor(targets$Sample_Group)
  design = model.matrix(~0 + groups)
  colnames(design) = sub("groups", "", colnames(design))
  aWeights <- arrayWeights(copy,design=design)
  
  
  #Permutations
  Zmatrix <- lapply(1:nPermutations, function(i) {
    print(i)
    set.seed(i)
    #shuffle data
    copy2 <- copy
    perm <- sample(1:length(targets$Sample_Group))
    aWeights2 <- aWeights[perm]
    targets2 <- targets[perm,]
    
    groups2 <- as.factor(targets2$Sample_Group)
    design2 =  model.matrix(~0 + groups2)
    colnames(design2) = sub("groups2", "", colnames(design2))
    fit_perm = lmFit(copy2, design2, weights=aWeights2)
    contrast.matrix2 = makeContrasts(res-sens,levels=design2)
    
    fit_perm <- contrasts.fit(fit_perm, contrast.matrix2)
    fit_perm <- eBayes(fit_perm)
    fit_perm$genes = array_exp$HGNC_ID
    sig_perm <- topTable(fit_perm, adjust.method="BH", coef=1,number = Inf,sort.by="none",resort.by=NULL)
    
    stat = data.frame(statistic=sig_perm$t,stringsAsFactors = F)
    rownames(stat)=sig_perm$ID
    return(stat)
    #get Z scores
    
    
  })
  Zmatrix <- do.call("cbind", Zmatrix)
  colnames(Zmatrix) <- 1:nPermutations
  return(Zmatrix)
}

nullmodel <- limmaNull( array_exp[,idy], targets, nPermutations=1000)
doxoSig = list(statistic = sig2_noSort$t, p.value=sig2_noSort$P.Value  )
names(doxoSig$statistic)=names(doxoSig$p.value)=sig2_noSort$ID
save(sig2_noSort,doxoSig,nullmodel,file = "../data/cellLineDoxoSignature.RData")

# Suppfig2
library(ggplot2)

exp_raw = array_exp[rownames(sig2)[which(sig2$P.Value<0.01)],]
exp_raw2 =apply(exp_raw[,-c(1)],1,as.numeric)
exp_imp = impute::impute.knn(exp_raw2)$data
rownames(exp_imp)= colnames(exp_raw)[-c(1)]



annotation_row = responses_raw
annotation_row$type = as.factor(ifelse(annotation_row$Doxo>qt[2],"sens",ifelse(annotation_row$Doxo<qt[1],"res",ifelse(is.na(annotation_row),NA,"int"))))
rownames(annotation_row)=annotation_row$validNames
idrem = which(is.na(annotation_row$Doxo))


tmp_exp = exp_imp[rownames(exp_imp) %in% annotation_row$validNames[-idrem],]
sort_doxo = annotation_row$validNames[annotation_row$validNames %in% rownames(tmp_exp)][order(annotation_row$Doxo[annotation_row$validNames %in% rownames(tmp_exp)])]
rem_1 = which(annotation_row$type[annotation_row$validNames %in% rownames(tmp_exp)][order(annotation_row$Doxo[annotation_row$validNames %in% rownames(tmp_exp)])] %in% c("res","sens"))
pdf("../plots/FigS3PanelA_heatmap_signature_dox.pdf",width = 7,height = 5)
pheatmap::pheatmap(tmp_exp[sort_doxo[rem_1],],scale = "column",annotation_row = annotation_row[,c(2,5),drop=F],show_colnames=F,
                   clustering_distance_rows = "correlation",clustering_distance_cols = "euclidean",clustering_method="average",
                   cluster_rows = F,color=colorRampPalette(c("navy", "white", "firebrick3"))(20) )
dev.off()

up <- sum(sig2$P.Value <= 0.01 & sig2$logFC >= 0.8)
down <- sum(sig2$P.Value <= 0.01 & sig2$logFC <= -0.8)
sig2$threshold <- as.factor(sig2$P.Value <= 0.05 & abs(sig2$logFC) >= 0.5) #for visualisation
sig2$logPval = -log10(sig2$P.Value)
library(ggrastr)
ggplot(data=sig2, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) + theme_bw()+scale_color_manual(values=c("navy", "firebrick3"),labels=c("Non significant","Sig"))+
  xlim(c(-3, 3)) + ylim(c(0, 7)) +
  xlab("log2 fold change") + ylab("-log10 p-value") 
ggsave("../plots/FigS3a_volcanodoxoSignature.pdf")
