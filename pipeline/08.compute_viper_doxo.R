# Compute VIPER enrichment

library(viper)
load("../data/chromatin_regulon.RData")
load("../data/complexes_v1.2_BRCA.RData")
load("../data/breast_doxo_resistance_data/saved_DES_aac.RData")

load("../data/cellLineDoxoSignature.RData")

mrs.a <- msviper(doxoSig[["statistic"]], chromatin_regulon, nullmodel, verbose =T)
p1 <- sort(mrs.a$es$p.value)
MRs <- names(p1)[p1 < 0.1]
MR_mc_heiser = MRs

MR_mc_heiser.p=mrs.a$es$p.value
MR_mc_heiser.nes=mrs.a$es$nes

MR_viper = function(signature, regulon){
  sig = signature$sig$t
  names(sig)=signature$sig$ID
  NM = as.matrix(signature$nullModel)
  rownames(NM)=signature$sig$ID
  mrs.a <- msviper(sig, regulon,NM , verbose =T)
  p1 <- sort(mrs.a$es$p.value)
  MRs <- names(p1)[p1 < 0.1]  
  list(MR=MRs,pval=mrs.a$es$p.value,nes=mrs.a$es$nes)
}





# CCLE_mc_CTRPv2
vipR = MR_viper(DE_CCLE_mc_CTRPv2_mc,chromatin_regulon)
MR_CCLE_mc_CTRPv2 = vipR$MR
MR_CCLE_mc_CTRPv2.p = vipR$pval
MR_CCLE_mc_CTRPv2.nes = vipR$nes

# CCLE_mc_finn_m
vipR = MR_viper(DE_CCLE_mc_finn_mc,chromatin_regulon)
MR_CCLE_mc_finn = vipR$MR
MR_CCLE_mc_finn.p = vipR$pval
MR_CCLE_mc_finn.nes = vipR$nes

# CCLE_mc_gray
vipR = MR_viper(DE_CCLE_mc_gray_mc,chromatin_regulon)
MR_CCLE_mc_gray = vipR$MR
MR_CCLE_mc_gray.p = vipR$pval
MR_CCLE_mc_gray.nes = vipR$nes

# CCLE_mc_GSDC
vipR = MR_viper(DE_CCLE_mc_GSDC1000_mc,chromatin_regulon)
MR_CCLE_mc_GSDC = vipR$MR
MR_CCLE_mc_GSDC.p = vipR$pval
MR_CCLE_mc_GSDC.nes =vipR$nes

# CCLE_rnasq_CTRP2
vipR= MR_viper(DE_CCLE_rnaseq_CTRPv2_rnaseq,chromatin_regulon)
MR_CCLE_rnasq_CTRP2 = vipR$MR
MR_CCLE_rnasq_CTRP2.p = vipR$pval
MR_CCLE_rnasq_CTRP2.nes = vipR$nes

# CCLE_rnasq_finn
vipR = MR_viper(DE_CCLE_rnaseq_finn_rnaseq,chromatin_regulon)
MR_CCLE_rnasq_finn = vipR$MR
MR_CCLE_rnasq_finn.p = vipR$pval
MR_CCLE_rnasq_finn.nes =  vipR$nes

# CCLE_rnasq_gray
vipR =  MR_viper(DE_CCLE_rnaseq_gray_rnaseq,chromatin_regulon)
MR_CCLE_rnasq_Gray = vipR$MR
MR_CCLE_rnasq_Gray.p = vipR$pval
MR_CCLE_rnasq_Gray.nes = vipR$nes


# rnaseeq_gcsi
vipR = MR_viper(DE_rnaseq_GCSI,chromatin_regulon)
MR_rnasq_gcsi = vipR$MR
MR_rnasq_gcsi.p = vipR$pval
MR_rnasq_gcsi.nes = vipR$nes

# rnaseq heiser
vipR = MR_viper(DE_heiser_rnaseq,chromatin_regulon)
MR_rnasq_heiser = vipR$MR
MR_rnasq_heiser.p = vipR$pval
MR_rnasq_heiser.nes = vipR$nes

combine.pval = data.frame(genes=names(chromatin_regulon),stringsAsFactors = F)
combine.pval$heiser_mc.p = MR_mc_heiser.p[match(combine.pval$genes,names(MR_mc_heiser.p))]
combine.pval$heiser_rnaseq.p = MR_rnasq_heiser.p[match(combine.pval$genes,names(MR_rnasq_heiser.p))]
combine.pval$MR_CCLE_mc_CTRPv2.p = MR_CCLE_mc_CTRPv2.p[match(combine.pval$genes,names(MR_CCLE_mc_CTRPv2.p))]
combine.pval$MR_CCLE_mc_finn.p = MR_CCLE_mc_finn.p[match(combine.pval$genes,names(MR_CCLE_mc_finn.p))]
combine.pval$MR_CCLE_mc_gray.p = MR_CCLE_mc_gray.p[match(combine.pval$genes,names(MR_CCLE_mc_gray.p))]
combine.pval$MR_CCLE_mc_GSDC.p = MR_CCLE_mc_GSDC.p[match(combine.pval$genes,names(MR_CCLE_mc_GSDC.p))]
combine.pval$MR_CCLE_rnasq_CTRP2.p = MR_CCLE_rnasq_CTRP2.p[match(combine.pval$genes,names(MR_CCLE_rnasq_CTRP2.p))]
combine.pval$MR_CCLE_rnasq_finn.p = MR_CCLE_rnasq_finn.p[match(combine.pval$genes,names(MR_CCLE_rnasq_finn.p))]
combine.pval$MR_CCLE_rnasq_Gray.p = MR_CCLE_rnasq_Gray.p[match(combine.pval$genes,names(MR_CCLE_rnasq_Gray.p))]
combine.pval$MR_rnasq_gcsi.p = MR_rnasq_gcsi.p[match(combine.pval$genes,names(MR_rnasq_gcsi.p))]

combine.nes = data.frame(genes=names(chromatin_regulon),stringsAsFactors = F)
combine.nes$heiser_mc.nes = MR_mc_heiser.nes[match(combine.nes$genes,names(MR_mc_heiser.nes))]
combine.nes$heiser_rnaseq.nes = MR_rnasq_heiser.nes[match(combine.nes$genes,names(MR_rnasq_heiser.nes))]
combine.nes$MR_CCLE_mc_CTRPv2.nes = MR_CCLE_mc_CTRPv2.nes[match(combine.nes$genes,names(MR_CCLE_mc_CTRPv2.nes))]
combine.nes$MR_CCLE_mc_finn.nes = MR_CCLE_mc_finn.nes[match(combine.nes$genes,names(MR_CCLE_mc_finn.nes))]
combine.nes$MR_CCLE_mc_gray.nes = MR_CCLE_mc_gray.nes[match(combine.nes$genes,names(MR_CCLE_mc_gray.nes))]
combine.nes$MR_CCLE_mc_GSDC.nes = MR_CCLE_mc_GSDC.nes[match(combine.nes$genes,names(MR_CCLE_mc_GSDC.nes))]
combine.nes$MR_CCLE_rnasq_CTRP2.nes = MR_CCLE_rnasq_CTRP2.nes[match(combine.nes$genes,names(MR_CCLE_rnasq_CTRP2.nes))]
combine.nes$MR_CCLE_rnasq_finn.nes = MR_CCLE_rnasq_finn.nes[match(combine.nes$genes,names(MR_CCLE_rnasq_finn.nes))]
combine.nes$MR_CCLE_rnasq_Gray.nes = MR_CCLE_rnasq_Gray.nes[match(combine.nes$genes,names(MR_CCLE_rnasq_Gray.nes))]
combine.nes$MR_rnasq_gcsi.nes = MR_rnasq_gcsi.nes[match(combine.nes$genes,names(MR_rnasq_gcsi.nes))]

nsamples = c(46,44,35,18,32,50,35,18,32,27)
combine.pval$combined =  apply(combine.pval[,-c(1)],1,function(x) survcomp::combine.test(x,nsamples,na.rm=T))
#combine.pval$combined1 =  apply(combine.pval[,-c(1)],1,function(x) survcomp::combine.test(x,na.rm=T))

print("Viper consensus enriched genes for doxorubicin signature")
print(combine.pval$genes[which(combine.pval$combined<0.05)])

save(MR_mc_heiser,MR_mc_heiser,MR_rnasq_heiser,MR_CCLE_mc_gray,MR_CCLE_mc_finn,MR_CCLE_mc_CTRPv2,MR_CCLE_mc_GSDC,MR_CCLE_rnasq_Gray,MR_CCLE_rnasq_finn,MR_CCLE_rnasq_CTRP2,MR_rnasq_gcsi,
     combine.pval,combine.nes,
     file="../data/cell_line_dataset_VIPER_datasets_aac.RData")


