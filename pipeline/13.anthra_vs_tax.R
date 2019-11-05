#13 surv anthra vs TAX with MKI67
# define which samples got anthra vs anthra+taxanes vs taxanes only vs CMF

load("../data/complexes_v1.2_BRCA.RData")
library(Biobase)
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}
# NEW DEFINITION
# Tax taxane only
# TaxAnt taxane + antra
# CMF 
# Ant antra only
# None no chemo


# load KAO
load("../data/clinical/rmas/GSE20685/rmaData.RData")
# 88 CAF vs 61 CMFs vs 119 others and 54 none
table(pData(eset_rma)$CX,pData(eset_rma)$neoadjuvant,useNA="a")

KAO = eset_rma
table(pData(KAO)$CX,pData(KAO)$anthra,useNA="a")
# anthraNA means noCT or no iidea
pData(KAO)$anthra[pData(KAO)$CX=="no"]=F
pData(KAO)$anthra[pData(KAO)$CX=="others"]=F

pData(KAO)$chemType = NA
pData(KAO)$chemType[which(pData(KAO)$CX == "others")]="Tax"
pData(KAO)$chemType[which(pData(KAO)$CX == "CAF")]="Ant"
pData(KAO)$chemType[which(pData(KAO)$CX == "CMF")]="CMF"
pData(KAO)$chemType[which(pData(KAO)$CX == "no")]="None"


# load INR
load("../data/clinical/rmas/GSE45255/rmaData.RData")
# antras, CMfs, tam only, tax only

INR = eset_rma
table(pData(INR)$CX,pData(INR)$anthra,useNA="a")
# false = none or other CN or tamox

pData(INR)$chemType = NA
pData(INR)$chemType[which(as.character(pData(INR)$CX) %in% c("treatment type: AC","treatment type: ACx4","treatment type: ACx4; Arimidex/Tamoxifen",
                                                             "treatment type: ACx4 cycles","treatment type: ACx4; Tamoxifen","treatment type: ACx6; Tamoxifen",
                                                             "treatment type: CAF","treatment type: CAFx6 cycles","treatment type: ECx4; Tamoxifen",
                                                             "treatment type: Tam+AC","treatment type: Tam+ACx4 cycles","treatment type: Tam+CAFx2 cycles",
                                                             "treatment type: Tam+CAFx6 cycles") )]="Ant"
pData(INR)$chemType[which(as.character(pData(INR)$CX) %in% c("treatment type: ACx8 cycles; taxol","treatment type: Tam+ACx8;taxol") )]="TaxAnt"
pData(INR)$chemType[which(as.character(pData(INR)$CX) %in% c("treatment type: anastrozole;taxol") )]="Tax"
pData(INR)$chemType[which(as.character(pData(INR)$CX) %in% c("treatment type: Arimidex+CMFx2 cycles","treatment type: CMF","treatment type: CMFx2; Tamoxifen",
                                                             "treatment type: CMFx6","treatment type: CMFx6; Arimidex/Tamoxifen","treatment type: CMFx6; Tamoxifen",
                                                             "treatment type: Tam+CMFx6 cycles") )]="CMF"
pData(INR)$chemType[which(as.character(pData(INR)$CX) %in%  c("treatment type: none","treatment type: Arimidex/Tamoxifen","treatment type: tamoxifen","treatment type: Tamoxifen",
                                                              "treatment type: tamoxifen;goserelin","treatment type: Tamoxifen+Zoledronic acid") )]="None"


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

pData(STK)$chemType = NA
pData(STK)$chemType[which(is.na(pData(STK)$CX))]="Ant"
pData(STK)$chemType[which(pData(STK)$CX == "CMF")]="CMF"
pData(STK)$chemType[which(pData(STK)$CX == "none")]="None"


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

pData(UPS)$chemType = "None"

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


pData(MAIRE)$chemType = NA
pData(MAIRE)$chemType[which(pData(MAIRE)$CX %in% c("Anthracyclines+Taxanes") )]="TaxAnt"
pData(MAIRE)$chemType[which(pData(MAIRE)$CX  %in% c( "Anthracyclines") )]="Ant"
pData(MAIRE)$chemType[which(pData(MAIRE)$CX == "Taxanes")]="Tax"

source("helper.R")


mergedES = mergeNONE(list(KAO,INR,STK,UPS,MAIRE))

pheno = pData(mergedES)
library(plyr)
pheno2 = colwise(function(x) plyr::revalue(x,c(" TRUE"="TRUE")))(pheno)
rownames(pheno2) = rownames(pheno)

pheno2$clinicSubtype = factor(ifelse(pheno2$her2=="TRUE","Her2+",
                                     ifelse(pheno2$er=="TRUE","ER+Her2-",
                                            ifelse(pheno2$pr=="FALSE","TNBC",NA))))



pheno =pheno2

batch = pData(mergedES)$cohort

merged_qn = mergedES
exprs(merged_qn)=preprocessCore::normalize.quantiles(exprs(mergedES))


mergedCOMBAT = removeBatchCombat(merged_qn,batch)
clinical_metacohort = mergedCOMBAT
# subtypes
expression_metacohort = t(exprs(mergedCOMBAT))
danot = mergedCOMBAT@featureData@data
danot$EntrezGene.ID =danot$ENTREZ_GENE_ID
danot$probe = rownames(danot)
danot$Gene.Symbol =danot$`Gene Symbol`

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



pheno_other = data.frame(rfs.t=pheno$os.t/365,rfs.e=pheno$os.e=="TRUE",chemo = pheno$chemType,age=pheno$age,
                         er = as.factor(pheno$er),pr = as.factor(pheno$pr),her2 = as.factor(pheno$her2),
                         lymphNodePos= pheno$lympNodePos,t.stage=factor(pheno$t_stage),cohort=pheno$cohort,
                         stringsAsFactors = F)

pheno_other$chemo[which(pheno_other$chemo %in% c("CMF","None","TaxAnt")) ]=NA

pheno_other_forTables = pheno_other[which(pheno_other$cohort %in% c("IRB/JNR/NUH","KAO","MAIRE","STK","UPP")),]
pheno_other_forTables$rfs.e =as.factor(pheno_other_forTables$rfs.e)


exp_raw = exprs(mergedCOMBAT)

pheno.df = na.omit(pheno_other)

pheno.df$cohort = droplevels(pheno.df$cohort)
pheno.df$chemo = droplevels(pheno.df$chemo)

pheno.df$MKI67_4 = rescale(exprs(mergedCOMBAT)[11409,rownames(pheno.df)])

eset_1_COMBAT2 = mergedCOMBAT[mergedCOMBAT@featureData@data$ENTREZ_GENE_ID %in% complexes$entrezgene,]

save(clinical_metacohort,pheno.df,file="../data/clinical_metacohort_anthra_vs_TAX.RData")

mulvarrfs.res = NULL
contrast.res = NULL
for(i in 1:nrow(eset_1_COMBAT2)){
  
  ezid = rownames(eset_1_COMBAT2)[i]
  
  X.b = cbind(pheno.df[,setdiff(colnames(pheno.df),"tmain")],exp_c=rescale(exprs(eset_1_COMBAT2)[i,rownames(pheno.df) ]))
  
  dd <- datadist(X.b)
  options(datadist='dd')
  
  try({
    
    cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*chemo+rcs(age,3)+er+pr+her2+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b)
    
    anv = anova(cxmodel2)
    pval = anv[14,3]
    coef = cxmodel2$coefficients
    min_lim = median(X.b$exp_c)-sd(X.b$exp_c)
    max_lim = median(X.b$exp_c)+sd(X.b$exp_c)
    
    min_lim = -1
    max_lim=1
    
    tmpCtrs = as.data.frame(rbind(contrast(cxmodel2,list(exp_c=min_lim,chemo="Ant"),list(exp_c=min_lim,chemo="Tax")),    
                                  contrast(cxmodel2,list(exp_c=max_lim,chemo="Ant"),list(exp_c=max_lim,chemo="Tax"))                                  
    ),stringAsFactors=F)
    
    tmpCtrs$gene = eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]
    tmpCtrs$probe = rownames(eset_1_COMBAT2)[i]
    contrast.res = rbind(contrast.res,tmpCtrs)
    
    smres = summary(cxmodel2)
    
    mulvarrfs.res = rbind(mulvarrfs.res,c(probe = ezid,entrez = as.character(eset_1_COMBAT2@featureData@data$ENTREZ_GENE_ID[i]), 
                                          gene=as.character(eset_1_COMBAT2@featureData@data$`Gene Symbol`[i]),
                                          pval.exp=anv[1,3],odds.exp=as.numeric(exp(coef)[1]),
                                          pval.anthr=as.numeric(anv[3,3]),odds.exp=as.numeric(exp(coef)[2]),
                                          pval.expAnth=as.numeric(anv[14,3]),odds.expAnth=as.numeric(exp(coef)[15])))
    
  })
  
  
  
  
}


mulvarrfs.res = data.frame(mulvarrfs.res,stringsAsFactors = F)
mulvarrfs.res$adjpval = p.adjust(as.numeric(mulvarrfs.res$pval.expAnth),method = "fdr")
mulvarrfs.res = mulvarrfs.res[order(as.numeric(mulvarrfs.res$pval.expAnth)),]


save(mulvarrfs.res,file="../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVStax_V2.RData")
save(contrast.res,file="../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVStax_V2.RData")

# table S6
mulvarrfs.res2 = mulvarrfs.res[mulvarrfs.res$entrez %in% complexes$entrezgene,]
write.table(mulvarrfs.res2,file="../anthra_vs_taxanes.tsv",quote = F,sep = "\t",row.names = F)


