# generate cell line signatures

breast_info = read.table("../data/breast_doxo_resistance_data/breast_doxorubicin_stat_table.csv",sep=",",header = T,stringsAsFactors = T)

library(PharmacoGx)
availablePSets()

breast_res = split(breast_info,breast_info$Dataset)

CCLE <- downloadPSet("CCLE")
gCSI <- downloadPSet("gCSI")
GDSC1000  <- downloadPSet("GDSC1000")
FIMM <- downloadPSet("FIMM")
CTRPv2 <- downloadPSet("CTRPv2")

save(CCLE,gCSI,GDSC1000,FIMM,CTRPv2,file="../data/breast_doxo_resistance_data/saved_ds.RData")

load("../data/breast_doxo_resistance_data/saved_ds.RData")

CCLE_mc_exp =  Biobase::exprs(CCLE@molecularProfiles[["rna"]])
colnames(CCLE_mc_exp)=CCLE@molecularProfiles$rna$cellid
CCLE_rnaseq_exp =  Biobase::exprs(CCLE@molecularProfiles[["rnaseq"]])
colnames(CCLE_rnaseq_exp)=CCLE@molecularProfiles$rnaseq$cellid
CCLE_cellInfo = CCLE@cell
CCLE_cellID = CCLE_cellInfo$cellid
gCSI_rnaseq_exp = Biobase::exprs(gCSI@molecularProfiles[["rnaseq"]])

GDSC1000_mc_exp =  Biobase::exprs(GDSC1000@molecularProfiles[["rna"]])
colnames(GDSC1000_mc_exp)=GDSC1000@molecularProfiles$rna$cellid

# CCLE_mc-gray 32
breast_info.Gray = breast_info[breast_info$Dataset=="GRAY",]
z = intersect(breast_info.Gray$Cell.Line,colnames(CCLE_mc_exp))
drug_CCLE_mc_gray_mc_exp = breast_info.Gray[match(z,breast_info.Gray$Cell.Line),]
exp_CCLE_mc_gray_mc_exp = CCLE_mc_exp[,z]

z_ccle_mc_grey = z

# CCLE_mc-FIM 18
breast_info.FIMM = breast_info[breast_info$Dataset=="FIMM",]
z = intersect(breast_info.FIMM$Cell.Line,colnames(CCLE_mc_exp))
drug_CCLE_mc_fimm_mc_exp = breast_info.FIMM[match(z,breast_info.FIMM$Cell.Line),]
exp_CCLE_mc_fimm_mc_exp = CCLE_mc_exp[,z]

z_ccle_mc_fim = z

# CCLE_mc-CTRPv2 35
breast_info.CTRPv2 = breast_info[breast_info$Dataset=="CTRPv2",]
z = intersect(breast_info.CTRPv2$Cell.Line,colnames(CCLE_mc_exp))
drug_CCLE_mc_CTRPv2_mc_exp = breast_info.CTRPv2[match(z,breast_info.CTRPv2$Cell.Line),]
exp_CCLE_mc_CTRPv2_mc_exp = CCLE_mc_exp[,z]

z_ccle_mc_CRTPV2 = z

# GDSC1000 (mc) 50
breast_info.GSDC1000 = breast_info[breast_info$Dataset=="GDSC1000",]
z = intersect(breast_info.GSDC1000$Cell.Line,colnames(GDSC1000_mc_exp))
drug_mc_GSDC1000_mc_exp = breast_info.GSDC1000[match(z,breast_info.GSDC1000$Cell.Line),]
exp_mc_GSDC1000_mc_exp = GDSC1000_mc_exp[,z]

z_mc_gdsc1k = z

# GCSI (rnaseq) 27
breast_info.GCSI = breast_info[breast_info$Dataset=="gCSI",]
breast_info.GCSI = breast_info.GCSI[breast_info.GCSI$Cell.Line!="HDQ-P1",]
z = intersect(breast_info.GCSI$Cell.Line,colnames(gCSI_rnaseq_exp))
drug_rnaseq_GCSI_rnaseq_exp = breast_info.GCSI[match(z,breast_info.GCSI$Cell.Line),]
exp_rnaseq_GCSI_rnaseq_exp = gCSI_rnaseq_exp[,z]

z_rnaseq_gsci = z

# CCLE_rnaseq-gray 32
breast_info.Gray = breast_info[breast_info$Dataset=="GRAY",]
z = intersect(breast_info.Gray$Cell.Line,colnames(CCLE_rnaseq_exp))
drug_CCLE_rnaseq_gray_rnaseq_exp = breast_info.Gray[match(z,breast_info.Gray$Cell.Line),]
exp_CCLE_rnaseq_gray_rnaseq_exp = CCLE_rnaseq_exp[,z]

z_CCLE_rnaseq_gray = z

# CCLE_rnaseq-FIM 18
breast_info.FIMM = breast_info[breast_info$Dataset=="FIMM",]
z = intersect(breast_info.FIMM$Cell.Line,colnames(CCLE_rnaseq_exp))
drug_CCLE_rnaseq_fimm_mc_rnaseq = breast_info.FIMM[match(z,breast_info.FIMM$Cell.Line),]
exp_CCLE_rnaseq_fimm_mc_rnaseq = CCLE_rnaseq_exp[,z]

z_CCLE_rnaseq_fimm = z

# CCLE_rnaseq-CTRPv2  35
breast_info.CTRPv2 = breast_info[breast_info$Dataset=="CTRPv2",]
z = intersect(breast_info.CTRPv2$Cell.Line,colnames(CCLE_rnaseq_exp))
drug_CCLE_rnaseq_CTRPv2_rnaseq_exp = breast_info.CTRPv2[match(z,breast_info.CTRPv2$Cell.Line),]
exp_CCLE_rnaseq_CTRPv2_rnaseq_exp = CCLE_rnaseq_exp[,z]

z_CCLE_rnaseq_CTRPv2 = z


# heiser 44
responses_raw = read.table("../data/heiser/procesedResponseData.tsv.csv",header = T,sep = "\t",stringsAsFactors = F)
colnames(responses_raw)=c("Cell.Line","GI50","score")
rnaseq_exp = read.table("../data/heiser/DREAM7_DrugSensitivity1_RNAseq_quantification.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T)
mc_exp = read.table("../data/heiser/DREAM7_DrugSensitivity1_GeneExpression.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T)
colnames(rnaseq_exp)[3:7]=c("184A1", "184B5",  "21MT1", "21NT",  "600MPE")
rnaseq_exp2 = rnaseq_exp[,3:46]
rownames(rnaseq_exp2)=rnaseq_exp$Ensembl_ID
z = intersect(responses_raw$Cell.Line,colnames(rnaseq_exp2))
drug_heiser_rnaseq = responses_raw[match(z,responses_raw$Cell.Line),]
exp_heiser_rnaseq = rnaseq_exp2[,z]

z_rnaseq_heiser = z

require(org.Hs.eg.db)
des_raw = select(org.Hs.eg.db, as.character(rnaseq_exp$Ensembl_ID), c("ENTREZID","GENENAME","SYMBOL","ENSEMBL"), "ENSEMBL")


heiser_featInfo = data.frame(Symbol=rnaseq_exp$HGNC_ID,EntrezGeneId=des_raw$ENTREZID[match(rnaseq_exp$Ensembl_ID,des_raw$ENSEMBL)],Ensembl =rnaseq_exp$Ensembl_ID )


all_z = Reduce(union,
               list(z_ccle_mc_CRTPV2,z_ccle_mc_fim,z_ccle_mc_grey,z_CCLE_rnaseq_CTRPv2,z_CCLE_rnaseq_fimm,z_CCLE_rnaseq_gray,z_mc_gdsc1k,z_rnaseq_gsci,z_rnaseq_heiser)
               )


runLimma = function(exp_data,drug_data,sens,featInfo,bootstrap=1000){
  require(limma)
  
  if(!is.numeric(drug_data[,sens])){
    drug_data[,sens] = as.numeric(as.character(drug_data[,sens]))
  }
  
  qt = quantile(drug_data[,sens],c(1/3,2/3),na.rm=T)
  
  group = as.factor(ifelse(drug_data[,sens]>qt[2],"sens",ifelse(drug_data[,sens]<qt[1],"res",NA)))
  design = model.matrix(~0+group[!is.na(group)])
  colnames(design) =c("res","sens")
  exp_data = exp_data[,colnames(exp_data)[!is.na(group)]]
  targets = data.frame(Sample_Name = colnames(exp_data),Sample_Group = group[!is.na(group)],stringsAsFactors = F)
  
  aWeights <- arrayWeights(exp_data,design=design)
  
  fit = lmFit(exp_data, design, weights=aWeights)
  contrast.matrix = makeContrasts(res-sens,levels=design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  std.error <- fit2$stdev.unscaled * fit2$sigma
  
  fit2$genes = featInfo$Symbol
  fit2$entrez = featInfo$EntrezGeneId
  sig2 <- topTable(fit2, adjust.method="BH", coef=1,number = Inf,sort.by="none",resort.by=NULL,confint = T) 
  # (cir-cil)/3.92
  sig2$std_error = std.error
  sig2$sigma = fit2$sigma
  
  limmaNull <- function(copy, targets, nPermutations=1000){
    #copy = array_exp[,idy]
    
    
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
      #fit_perm$genes = array_exp$HGNC_ID
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
  
  nullmodel <- limmaNull( exp_data, targets, nPermutations=bootstrap)
  
  res = list()
  res$sig = sig2
  res$nullModel = nullmodel
  res$featInfor = featInfo
  return(res)  
  
}


runLimmaRNAseq = function(exp_data,drug_data,sens,featInfo,bootstrap=1000,norm=F){
  require(limma)
  
  if(!is.numeric(drug_data[,sens])){
    drug_data[,sens] = as.numeric(as.character(drug_data[,sens]))
  }
  
  qt = quantile(drug_data[,sens],c(1/3,2/3),na.rm=T)
  
  group = as.factor(ifelse(drug_data[,sens]>qt[2],"sens",ifelse(drug_data[,sens]<qt[1],"res",NA)))
  design = model.matrix(~0+group[!is.na(group)])
  colnames(design) =c("res","sens")
  exp_data = exp_data[,colnames(exp_data)[!is.na(group)]]
  if(norm){
    exp_data = voom(exp_data,design) 
  }
  
  targets = data.frame(Sample_Name = colnames(exp_data),Sample_Group = group[!is.na(group)],stringsAsFactors = F)
  
  #aWeights <- arrayWeights(exp_data,design=design)
  
  #fit = lmFit(exp_data, design, weights=aWeights)
  fit = lmFit(exp_data, design)
  contrast.matrix = makeContrasts(res-sens,levels=design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  fit2$genes = featInfo$Symbol
  fit2$entrez = featInfo$EntrezGeneId
  sig2 <- topTable(fit2, adjust.method="BH", coef=1,number = Inf,sort.by="none",resort.by=NULL,confint = T) 
  
  limmaNull <- function(copy, targets, nPermutations=1000){
    #copy = array_exp[,idy]
    
    
    #Estimates relative quality weights for each array in a multi-array experiment.
    groups <- as.factor(targets$Sample_Group)
    design = model.matrix(~0 + groups)
    colnames(design) = sub("groups", "", colnames(design))
    #aWeights <- arrayWeights(copy,design=design)
    
    
    #Permutations
    Zmatrix <- lapply(1:nPermutations, function(i) {
      print(i)
      set.seed(i)
      #shuffle data
      copy2 <- copy
      perm <- sample(1:length(targets$Sample_Group))
      #aWeights2 <- aWeights[perm]
      targets2 <- targets[perm,]
      
      groups2 <- as.factor(targets2$Sample_Group)
      design2 =  model.matrix(~0 + groups2)
      colnames(design2) = sub("groups2", "", colnames(design2))
      #fit_perm = lmFit(copy2, design2, weights=aWeights2)
      fit_perm = lmFit(copy2, design2)
      contrast.matrix2 = makeContrasts(res-sens,levels=design2)
      
      fit_perm <- contrasts.fit(fit_perm, contrast.matrix2)
      fit_perm <- eBayes(fit_perm)
      #fit_perm$genes = array_exp$HGNC_ID
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
  
  nullmodel <- limmaNull( exp_data, targets, nPermutations=bootstrap)
  
  res = list()
  res$sig = sig2
  res$nullModel = nullmodel
  res$featInfor = featInfo
  return(res)  
  
}

DE_CCLE_mc_gray_mc = runLimma(exp_CCLE_mc_gray_mc_exp,drug_CCLE_mc_gray_mc_exp,sens="AAC....",featInfo =CCLE@molecularProfiles$rna@featureData@data,bootstrap = 1000 )
DE_CCLE_mc_finn_mc = runLimma(exp_CCLE_mc_fimm_mc_exp,drug_CCLE_mc_fimm_mc_exp,sens="AAC....",featInfo =CCLE@molecularProfiles$rna@featureData@data ,bootstrap = 1000)
DE_CCLE_mc_CTRPv2_mc = runLimma(exp_CCLE_mc_CTRPv2_mc_exp,drug_CCLE_mc_CTRPv2_mc_exp,sens="AAC....",featInfo =CCLE@molecularProfiles$rna@featureData@data ,bootstrap = 1000)
DE_CCLE_mc_GSDC1000_mc = runLimma(exp_mc_GSDC1000_mc_exp,drug_mc_GSDC1000_mc_exp,sens="AAC....",featInfo =GDSC1000@molecularProfiles$rna@featureData@data ,bootstrap = 1000)
DE_rnaseq_GCSI = runLimmaRNAseq(exp_rnaseq_GCSI_rnaseq_exp,drug_rnaseq_GCSI_rnaseq_exp,sens = "AAC....",featInfo = gCSI@molecularProfiles$rnaseq@featureData@data,bootstrap = 1000)
DE_CCLE_rnaseq_gray_rnaseq = runLimmaRNAseq(exp_CCLE_rnaseq_gray_rnaseq_exp,drug_CCLE_rnaseq_gray_rnaseq_exp,sens="AAC....",featInfo =CCLE@molecularProfiles$rnaseq@featureData@data,bootstrap = 1000)
DE_CCLE_rnaseq_finn_rnaseq = runLimmaRNAseq(exp_CCLE_rnaseq_fimm_mc_rnaseq,drug_CCLE_rnaseq_fimm_mc_rnaseq,sens="AAC....",featInfo =CCLE@molecularProfiles$rnaseq@featureData@data ,bootstrap = 1000)
DE_CCLE_rnaseq_CTRPv2_rnaseq = runLimmaRNAseq(exp_CCLE_rnaseq_CTRPv2_rnaseq_exp,drug_CCLE_rnaseq_CTRPv2_rnaseq_exp,sens="AAC....",featInfo =CCLE@molecularProfiles$rnaseq@featureData@data ,bootstrap = 1000)
DE_heiser_rnaseq = runLimmaRNAseq(exp_heiser_rnaseq,drug_heiser_rnaseq,sens="GI50",featInfo =heiser_featInfo,norm = T,bootstrap = 1000)

DE_heiser_rnaseq_sve = DE_heiser_rnaseq
DE_heiser_rnaseq = DE_heiser_rnaseq_sve
save(DE_CCLE_mc_gray_mc,DE_CCLE_mc_finn_mc,DE_CCLE_mc_CTRPv2_mc,DE_CCLE_mc_CTRPv2_mc,
     DE_CCLE_mc_GSDC1000_mc,DE_rnaseq_GCSI,DE_CCLE_rnaseq_gray_rnaseq,DE_CCLE_rnaseq_finn_rnaseq,
     DE_CCLE_rnaseq_CTRPv2_rnaseq,DE_heiser_rnaseq,file="../data/breast_doxo_resistance_data/saved_DES_aac.RData")

