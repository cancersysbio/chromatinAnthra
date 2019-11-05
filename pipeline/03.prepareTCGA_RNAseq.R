library(DESeq2)
library(org.Hs.eg.db)



load('../data/02.GeneExp_gene_counts.RData')

exp_mat = MultiAssayExperiment::assay(data.exp)
rn = rownames(exp_mat)

Des = data.exp@rowRanges
coldata = data.exp@colData
dds <- DESeqDataSetFromMatrix(countData = exp_mat,
                              colData = coldata,
                              design= ~  shortLetterCode)

vsd <- vst(dds, blind=T)
exp.vst = MultiAssayExperiment::assay(vsd)
ts = substr(colnames(exp.vst),14,15)

exp.T.vst = exp.vst[,ts=="01"]
exp.N.vst = exp.vst[,which(ts=="11")]

tmpSamples = substr(colnames(exp.T.vst),1,15)
idDup = duplicated(tmpSamples)

exp.T.vst = exp.T.vst[,which(!idDup)]



save(exp.T.vst,exp.N.vst,Des,file="../data/vst.RData")

write.table(cbind(gene= rownames(exp.T.vst),exp.T.vst),file="../data/aracne_APtumor.exp",sep="\t",row.names = F,col.names = T,quote=F)

write.table(rownames(exp.T.vst),file="../data/aracne_APgeneList.exp",sep="\t",row.names = F,col.names = F,quote=F)

