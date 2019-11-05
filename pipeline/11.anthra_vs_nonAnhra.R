## anthra vs non anthra

load("../data/complexes_v1.2_BRCA.RData")
library(Biobase)
library(viridis)
library(ggrastr)

rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}


# load KAO
load("../data/clinical/rmas/GSE20685/rmaData.RData")
# 88 CAF vs 61 CMFs vs 119 others and 54 none
table(pData(eset_rma)$CX,pData(eset_rma)$neoadjuvant,useNA="a")

KAO = eset_rma
table(pData(KAO)$CX,pData(KAO)$anthra,useNA="a")
# anthraNA means noCT or no iidea
pData(KAO)$anthra[pData(KAO)$CX=="no"]=F
pData(KAO)$anthra[pData(KAO)$CX=="others"]=F


# load INR
load("../data/clinical/rmas/GSE45255/rmaData.RData")
# antras, CMfs, tam only, tax only

INR = eset_rma
table(pData(INR)$CX,pData(INR)$anthra,useNA="a")
# false = none or other CN or tamox

# load STK
load("../data/clinical/rmas/GSE1456/rmaData_96.RData")
STK_96 = eset_rma
STK_96 = STK_96[,1:159]
colnames(STK_96) = rownames(pData(STK_96))

load("../data/clinical/rmas/GSE1456/rmaData_97.RData")

STK_97 = eset_rma
STK_97 = STK_97[,1:159]
colnames(STK_97) = rownames(pData(STK_97))

z = intersect(rownames(STK_96),rownames(STK_97))
z97 = setdiff(rownames(STK_97),rownames(STK_96))
combined_exp = rbind(exprs(STK_96),exprs(STK_97)[z97,])
combined_fs = rbind(featureData(STK_96)@data,featureData(STK_97)@data[z97,])


STK = ExpressionSet(assayData = combined_exp)
featureData(STK)=AnnotatedDataFrame(combined_fs)
colnames(STK)= colnames(STK_97)
rownames(pData(STK_97))
pData(STK) = pData(STK_96)[1:159,]
STK@annotation="GPL97_GPL96"

# 90 CMF o 30 none 47 NA (likely CMF to)
table(pData(STK)$CX,pData(STK)$anthra,useNA="a")
# here false is CM or none

# load UPS
load("../data/clinical/rmas/GSE3494/rmaData_96.RData")
UPS_96 = eset_rma
colnames(UPS_96) = rownames(pData(UPS_96))


load("../data/clinical/rmas/GSE3494/rmaData_97.RData")
UPS_97 = eset_rma


z = intersect(rownames(UPS_96),rownames(UPS_97))
z97 = setdiff(rownames(UPS_97),rownames(UPS_96))
combined_exp = rbind(exprs(UPS_96),exprs(UPS_97)[z97,])
combined_fs = rbind(featureData(UPS_96)@data,featureData(UPS_97)@data[z97,])


UPS = ExpressionSet(assayData = combined_exp)
featureData(UPS)=AnnotatedDataFrame(combined_fs)
pData(UPS) = pData(UPS_96)
UPS@annotation="GPL97_GPL96"
table(pData(UPS)$CX,pData(UPS)$anthra,useNA="a")
# false no treatment


# load MAIRE
load("../data/clinical/rmas/GSE65194/rmaData.RData")
# 70+23 anthra, 30 no treatment
MAIRE = eset_rma
MAIRE = MAIRE[,rownames(pData(MAIRE))]

pData(MAIRE)$her2 = plyr::revalue(pData(MAIRE)$her2,c("no"="FALSE","yes"="TRUE"))

table(pData(MAIRE)$chemo)
levels(pData(MAIRE)$chemo)

table(pData(MAIRE)$CX,pData(MAIRE)$anthra,useNA="a")
# false NA or other or taxa
pData(MAIRE)$anthra[is.na(pData(MAIRE)$CX)]=NA

source("helper.R")


mergedES = mergeNONE(list(KAO,INR,STK,UPS,MAIRE))

pheno = pData(mergedES)
library(plyr)

pheno2 = colwise(function(x) plyr::revalue(x,c(" TRUE"="TRUE")))(pheno)
rownames(pheno2) = rownames(pheno)

pheno =pheno2

batch = pData(mergedES)$cohort

merged_qn = mergedES
exprs(merged_qn)=preprocessCore::normalize.quantiles(exprs(mergedES))


mergedCOMBAT = removeBatchCombat(merged_qn,batch)

clinical_metacohort = mergedCOMBAT

pheno = pData(mergedCOMBAT)

pheno2 = colwise(function(x) plyr::revalue(x,c(" TRUE"="TRUE")))(pheno)
rownames(pheno2) = rownames(pheno)

pheno2$clinicSubtype = factor(ifelse(pheno2$her2=="TRUE","Her2+",
                                    ifelse(pheno2$er=="TRUE","ER+Her2-",
                                           ifelse(pheno2$pr=="FALSE","TNBC",NA))))

pheno=pheno2
pheno$age = as.numeric(as.character(pheno$age))
pheno$lymphNodeNum = as.numeric(as.character(pheno$lymphNodeNum))
pheno$os.t = as.numeric(as.character(pheno$os.t))
pheno$rdfs.t = as.numeric(as.character(pheno$rdfs.t))
pheno$rfs.t = as.numeric(as.character(pheno$rfs.t))
pheno$size = as.numeric(as.character(pheno$size))
pheno$t_stage = factor(as.numeric(as.character(pheno$t_stage)))

pheno_other = data.frame(rfs.t=pheno$os.t/365,rfs.e=pheno$os.e=="TRUE",anthra=pheno$anthra,age=pheno$age,
                         er = as.factor(pheno$er),pr = as.factor(pheno$pr),her2 = as.factor(pheno$her2),
                         lymphNodePos= pheno$lympNodePos,t.stage=factor(pheno$t_stage),cohort=pheno$cohort,
                         stringsAsFactors = F)

pheno.df = na.omit(pheno_other)
pheno.df$cohort = droplevels(pheno.df$cohort)

library(rms)
library(survival)


pheno.df$MKI67_4 = rescale(exprs(mergedCOMBAT)[11409,rownames(pheno.df)])

eset_1_COMBAT2 = mergedCOMBAT[mergedCOMBAT@featureData@data$ENTREZ_GENE_ID %in% complexes$entrezgene,]

save(clinical_metacohort,pheno.df,file="../data/clinical_metacohort_anthra_vs_nonAnthra.RData")

mulvarrfs.res =mulvarrfs.res.erP=mulvarrfs.res.her2P=mulvarrfs.res.TNBC= NULL
contrast.res = contrast.res.erP=contrast.res.her2P=contrast.res.TNBC= NULL
for(i in 1:nrow(eset_1_COMBAT2)){
  
  ezid = rownames(eset_1_COMBAT2)[i]
  esr1 = rescale(exprs(eset_1_COMBAT2)[188,rownames(pheno.df) ])
  X.b = cbind(pheno.df[,setdiff(colnames(pheno.df),"tmain")],exp_c=rescale(exprs(eset_1_COMBAT2)[i,rownames(pheno.df) ]),esr1=esr1)
  
  dd <- datadist(X.b)
  options(datadist='dd')
  
  try({
    
    cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+er+pr+her2+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b)
    
    
    anv = anova(cxmodel2)
    pval = anv[14,3]
    coef = cxmodel2$coefficients
    min_lim = median(X.b$exp_c)-sd(X.b$exp_c)
    max_lim = median(X.b$exp_c)+sd(X.b$exp_c)
    
    min_lim = -1
    max_lim=1
    
    
    tmpCtrs = as.data.frame(rbind(contrast(cxmodel2,list(exp_c=min_lim,anthra=T),list(exp_c=min_lim,anthra=F)),
                                  contrast(cxmodel2,list(exp_c=max_lim,anthra=T),list(exp_c=max_lim,anthra=F))
    ),stringAsFactors=F)
    
    tmpCtrs$gene = eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]
    tmpCtrs$probe = rownames(eset_1_COMBAT2)[i]
    contrast.res = rbind(contrast.res,tmpCtrs)
    
    smres = summary(cxmodel2)
    
    mulvarrfs.res = rbind(mulvarrfs.res,c(probe = ezid,entrez = as.character(eset_1_COMBAT2@featureData@data$ENTREZ_GENE_ID[i]), 
                                          gene=as.character(eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]),
                                          pval.exp=anv[1,3],odds.exp=as.numeric(exp(coef)[1]),
                                          pval.anthr=as.numeric(anv[3,3]),odds.exp=as.numeric(exp(coef)[2]),
                                          pval.expAnth=as.numeric(anv[14,3]),odds.expAnth=as.numeric(exp(coef)[16])))
    
  })
  
  
  
  
  # ER+
  try({
    X.b.erP = X.b[which(X.b$er=="TRUE" & X.b$her2=="FALSE"),]
    
    cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+pr+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b.erP)
    
    
    anv = anova(cxmodel2)
    pval = anv[12,3]
    coef = cxmodel2$coefficients
    min_lim = median(X.b.erP$exp_c)-sd(X.b.erP$exp_c)
    max_lim = median(X.b.erP$exp_c)+sd(X.b.erP$exp_c)
    
    min_lim = -1
    max_lim=1
    
    
    tmpCtrs = as.data.frame(rbind(contrast(cxmodel2,list(exp_c=min_lim,anthra=T),list(exp_c=min_lim,anthra=F)),
                                  contrast(cxmodel2,list(exp_c=max_lim,anthra=T),list(exp_c=max_lim,anthra=F))
    ),stringAsFactors=F)
    
    
    tmpCtrs$gene = eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]
    tmpCtrs$probe = rownames(eset_1_COMBAT2)[i]
    contrast.res.erP = rbind(contrast.res.erP,tmpCtrs)
    
    smres = summary(cxmodel2)
    
    mulvarrfs.res.erP = rbind(mulvarrfs.res.erP,c(probe = ezid,entrez = as.character(eset_1_COMBAT2@featureData@data$ENTREZ_GENE_ID[i]), 
                                                  gene=as.character(eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]),
                                                  pval.exp=anv[1,3],odds.exp=as.numeric(exp(coef)[1]),
                                                  pval.anthr=as.numeric(anv[3,3]),odds.exp=as.numeric(exp(coef)[2]),
                                                  pval.expAnth=as.numeric(anv[12,3]),odds.expAnth=as.numeric(exp(coef)[14])))
    
  })
  
  
  # Her2+
  try({
    X.b.her2P = X.b[which(X.b$her2=="TRUE"),]
    cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+er+pr+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b.her2P)
        
    anv = anova(cxmodel2)    
    coef = cxmodel2$coefficients
    min_lim = median(X.b.her2P$exp_c)-sd(X.b.her2P$exp_c)
    max_lim = median(X.b.her2P$exp_c)+sd(X.b.her2P$exp_c)
    
    min_lim = -1
    max_lim=1
        
    tmpCtrs = as.data.frame(rbind(contrast(cxmodel2,list(exp_c=min_lim,anthra=T),list(exp_c=min_lim,anthra=F)),
                                  contrast(cxmodel2,list(exp_c=max_lim,anthra=T),list(exp_c=max_lim,anthra=F))
    ),stringAsFactors=F)
    
    
    tmpCtrs$gene = eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]
    tmpCtrs$probe = rownames(eset_1_COMBAT2)[i]
    contrast.res.her2P = rbind(contrast.res.her2P,tmpCtrs)
    
    smres = summary(cxmodel2)
    
    mulvarrfs.res.her2P = rbind(mulvarrfs.res.her2P,c(probe = ezid,entrez = as.character(eset_1_COMBAT2@featureData@data$ENTREZ_GENE_ID[i]), 
                                                      gene=as.character(eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]),
                                                      pval.exp=anv[1,3],odds.exp=as.numeric(exp(coef)[1]),
                                                      pval.anthr=as.numeric(anv[3,3]),odds.exp=as.numeric(exp(coef)[2]),
                                                      pval.expAnth=as.numeric(anv[13,3]),odds.expAnth=as.numeric(exp(coef)[15])))
    
  })
  
  # TNBC
  try({
    X.b.TNBC =X.b[which(X.b$her2=="FALSE"&X.b$er=="FALSE"&X.b$pr=="FALSE"),] 
    cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b.TNBC)
    
    anv = anova(cxmodel2)    
    coef = cxmodel2$coefficients
    min_lim = median(X.b.TNBC$exp_c)-sd(X.b.TNBC$exp_c)
    max_lim = median(X.b.TNBC$exp_c)+sd(X.b.TNBC$exp_c)
    
    min_lim = -1
    max_lim=1
    
    
    tmpCtrs = as.data.frame(rbind(contrast(cxmodel2,list(exp_c=min_lim,anthra=T),list(exp_c=min_lim,anthra=F)),
                                  contrast(cxmodel2,list(exp_c=max_lim,anthra=T),list(exp_c=max_lim,anthra=F))
    ),stringAsFactors=F)
    
    
    tmpCtrs$gene = eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]
    tmpCtrs$probe = rownames(eset_1_COMBAT2)[i]
    contrast.res.TNBC = rbind(contrast.res.TNBC,tmpCtrs)
    
    smres = summary(cxmodel2)
    
    mulvarrfs.res.TNBC = rbind(mulvarrfs.res.TNBC,c(probe = ezid,entrez = as.character(eset_1_COMBAT2@featureData@data$ENTREZ_GENE_ID[i]), 
                                                    gene=as.character(eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]),
                                                    pval.exp=anv[1,3],odds.exp=as.numeric(exp(coef)[1]),
                                                    pval.anthr=as.numeric(anv[3,3]),odds.exp=as.numeric(exp(coef)[2]),
                                                    pval.expAnth=as.numeric(anv[11,3]),odds.expAnth=as.numeric(exp(coef)[13])))
    
  })
  
  
}


mulvarrfs.res = data.frame(mulvarrfs.res,stringsAsFactors = F)
mulvarrfs.res$adjpval = p.adjust(as.numeric(mulvarrfs.res$pval.expAnth),method = "fdr")
mulvarrfs.res = mulvarrfs.res[order(as.numeric(mulvarrfs.res$pval.expAnth)),]


mulvarrfs.res.erP = data.frame(mulvarrfs.res.erP,stringsAsFactors = F)
mulvarrfs.res.erP$adjpval = p.adjust(as.numeric(mulvarrfs.res.erP$pval.expAnth),method = "fdr")
mulvarrfs.res.erP = mulvarrfs.res.erP[order(as.numeric(mulvarrfs.res.erP$pval.expAnth)),]

mulvarrfs.res.her2P = data.frame(mulvarrfs.res.her2P,stringsAsFactors = F)
mulvarrfs.res.her2P$adjpval = p.adjust(as.numeric(mulvarrfs.res.her2P$pval.expAnth),method = "fdr")
mulvarrfs.res.her2P = mulvarrfs.res.her2P[order(as.numeric(mulvarrfs.res.her2P$pval.expAnth)),]

mulvarrfs.res.TNBC= data.frame(mulvarrfs.res.TNBC,stringsAsFactors = F)
mulvarrfs.res.TNBC$adjpval = p.adjust(as.numeric(mulvarrfs.res.TNBC$pval.expAnth),method = "fdr")
mulvarrfs.res.TNBC= mulvarrfs.res.TNBC[order(as.numeric(mulvarrfs.res.TNBC$pval.expAnth)),]



save(mulvarrfs.res,
     mulvarrfs.res.erP,mulvarrfs.res.her2P,mulvarrfs.res.TNBC,file="../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")
save(contrast.res,
     contrast.res.erP,contrast.res.her2P,contrast.res.TNBC,file="../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")

# supp table S4
mulvarrfs.res2 = mulvarrfs.res[mulvarrfs.res$entrez %in% complexes$entrezgene,]
write.table(mulvarrfs.res2,file="../anthra_vs_nonAnthra.tsv",quote = F,sep = "\t",row.names = F)

# supp table S7
mulvarrfs.res.erP2 = mulvarrfs.res.erP[mulvarrfs.res.erP$entrez %in% complexes$entrezgene,]
write.table(mulvarrfs.res.erP2,file="../anthra_vs_nonAnthra_ERPos.tsv",quote = F,sep = "\t",row.names = F)

# supp table S8
mulvarrfs.res.her2P2 = mulvarrfs.res.her2P[mulvarrfs.res.her2P$entrez %in% complexes$entrezgene,]
write.table(mulvarrfs.res.her2P2,file="../anthra_vs_nonAnthra_Her2Pos.tsv",quote = F,sep = "\t",row.names = F)

# supp Table S9
mulvarrfs.res.TNBC2 = mulvarrfs.res.TNBC[mulvarrfs.res.TNBC$entrez %in% complexes$entrezgene,]
write.table(mulvarrfs.res.TNBC2,file="../anthra_vs_nonAnthra_TNBC.tsv",quote = F,sep = "\t",row.names = F)



