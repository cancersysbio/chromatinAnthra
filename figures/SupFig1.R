# Supp figure 1

library(maftools)
load("/data2/TCGA_updated/PANCAN_33/mut/mutect2_mut_amp_del_MAF.RData")
load("complexes_v1.2_BRCA.RData")



pancan_MAF =subsetMaf(pancan_MAF_amp_del,query = "!Variant_Classification %in% c('Amp','Del')",mafObj = T,includeSyn = F)
mcm1 = mutCountMatrix(pancan_MAF,removeNonMutated = F)>0

mcm1 = mcm1[rownames(mcm1) %in% complexes$hgnc_symbol,]
BRCAness = cs1$subtype[match(colnames(mcm1),cs1$Tumor_Sample_Barcode)]=="TCGA-BRCA"

mcm2 = data.frame(t(mcm1),BRCAness,stringsAsFactors = F)
mcm3 = ddply(mcm2,.(BRCAness),function(x) apply(x,2,function(y) sum(y,na.rm=T)),.progress = "text")


mcm4 = data.frame(gene= colnames(mcm3),BRCA = as.numeric(mcm3[2,,drop=T]),PANCAN=as.numeric(mcm3[1,,drop=T]),stringsAsFactors = F)
mcm4 = mcm4[mcm4$gene!="BRCAness",]
mcm4$BRCA = mcm4$BRCA/985
mcm4$PANCAN = mcm4$PANCAN/9097
mcm4$subtype = complexes$name[match(mcm4$gene,complexes$hgnc_symbol)]
mcm4$color = complexes$color[match(mcm4$gene,complexes$hgnc_symbol)]

mcm4$chromatin_gene = ifelse(mcm4$subtype=="BRCA","Breast","Chromatin")

gois_mut = c("KMT2C","KMT2D","ARID1A","ATRX","ATM","ARID1B","FOXA1","KMT2A","RB1","ASH1L","KMT2B","ARID4B","BAZ2B","SMARCA4",
             "HIST1H3B","SMYD3","BRCA1","KMT2E","CTCF","SRCAP")

ggplot(subset(mcm4,chromatin_gene=="Chromatin"),aes(x=BRCA,y=PANCAN,label=gene))+
  geom_point(data=subset(mcm4,chromatin_gene=="Chromatin" & !gene %in% gois_mut),aes(x=BRCA,y=PANCAN,label=gene,alpha=0.5) )+
  geom_text(data = subset(mcm4,gene %in% gois_mut),aes(x=BRCA,y=PANCAN,label=gene,col="#F8766D"))+
  geom_abline(intercept = 0,slope=1)+theme_bw()+guides(alpha=FALSE,color=F)+
  scale_y_continuous(labels = scales::percent)+ scale_x_continuous(labels = scales::percent)