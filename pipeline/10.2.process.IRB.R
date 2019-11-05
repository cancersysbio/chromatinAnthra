library(GEOquery)

rawsetL = getGEO("GSE45255")
save(rawsetL,file = "./data/clinical//GSE45255_raw.RData")




expOrig = exprs(rawsetL[[1]])
prExp = expOrig["208305_at",]
her2Exp = expOrig["216836_s_at",]
erExp = expOrig["205225_at",]

getThreshold = function(exp){
  intersect <- function(m1, s1, m2, s2, prop1, prop2){
    
    B <- (m1/s1^2 - m2/s2^2)
    A <- 0.5*(1/s2^2 - 1/s1^2)
    C <- 0.5*(m2^2/s2^2 - m1^2/s1^2) - log((s1/s2)*(prop2/prop1))
    
    (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
  }
  library(flexmix)
  kk = flexmix(exp~1,k=2,model = list(FLXMRglm(exp ~ ., family = "gaussian"), FLXMRglm(exp ~ ., family = "gaussian")))
  mu1 = parameters(kk)[[1]][1,1]
  mu2 = parameters(kk)[[1]][1,2]
  sigma1 = parameters(kk)[[1]][2,1]
  sigma2 = parameters(kk)[[1]][2,2]
  hist(exp,breaks = 100,col = kk@cluster)
  intersect(mu1,sigma1,mu2,sigma2,table(clusters(kk))[1]/length(exp),table(clusters(kk))[2]/length(exp))[2]
  
}

getThreshold(prExp)
thresholdPR = 5.12601
getThreshold(her2Exp)
thresholdher2 = 11.15137
getThreshold(erExp)
thresholder =11.34597
pr_probe = prExp>thresholdPR
her2_probe = her2Exp>thresholdher2
er_probe = erExp>thresholder

phenoData = pData(rawsetL[[1]])

covariate.df= data.frame(title=phenoData$title, stringsAsFactors = F)
rownames(covariate.df)=rownames(phenoData)
covariate.df$adjuvant = T
covariate.df$adjuvant[phenoData$characteristics_ch1.8=="adjuvant treated? (0=no, 1=yes): 0"]=F
covariate.df$adjuvant[phenoData$characteristics_ch1.8=="adjuvant treated? (0=no, 1=yes): NA"]=NA

covariate.df$neoadjuvant = F

covariate.df$adjNONE = F
covariate.df$adjNONE[phenoData$characteristics_ch1.8=="adjuvant treated? (0=no, 1=yes): 0"] = T
covariate.df$adjNONE[phenoData$characteristics_ch1.8=="adjuvant treated? (0=no, 1=yes): NA"] = NA

covariate.df$adjNA = F
covariate.df$adjNA[phenoData$characteristics_ch1.8=="adjuvant treated? (0=no, 1=yes): NA"] = T


covariate.df$geoAc = phenoData$geo_accession
covariate.df$cohort = "IRB/JNR/NUH"

covariate.df$age =as.numeric(gsub("[^0-9]", "", phenoData$characteristics_ch1.6))
ER_1 = phenoData$characteristics_ch1.1=="er status: ER+"
table(ER,ER_1)
table(er_probe,ER_1)

PR_1 = phenoData$characteristics_ch1.2=="pgr status: PgR+"
table(PR,PR_1)
table(pr_probe,PR_1)

Her2_1 = phenoData$characteristics_ch1.3=="her2 status: He+"
table(HER2,Her2_1)
table(her2_probe,Her2_1)

covariate.df$er = ER_1
covariate.df$pr = PR_1
covariate.df$her2 = Her2_1
covariate.df$hist = phenoData$characteristics_ch1.7
grade = rep(NA,nrow(covariate.df))
grade[phenoData$characteristics_ch1.4=="histological grade: G1"] = 1
grade[phenoData$characteristics_ch1.4=="histological grade: G2"] = 2
grade[phenoData$characteristics_ch1.4=="histological grade: G3"] = 3

covariate.df$grade = grade


lymphNodeNum = NA
size= as.numeric(gsub("[^0-9]", "", phenoData$characteristics_ch1.5))


covariate.df$stage = NA
covariate.df$t_stage = ifelse(size<20,1,ifelse(size<50,2,3))
covariate.df$m_stage = NA
covariate.df$n_stage= NA
covariate.df$lympNodePos = phenoData$characteristics_ch1=="ln status: LN+"
covariate.df$lymphNodeNum = lymphNodeNum
covariate.df$size= size
covariate.df$CX = phenoData$characteristics_ch1.11
covariate.df$RX= NA
covariate.df$HX = phenoData$characteristics_ch1.9
covariate.df$anthra = F
idAntra = which(phenoData$characteristics_ch1.11   %in% c("treatment type: AC","treatment type: ACx3/CMFx1","treatment type: ACx4",
                                                         "treatment type: ACx4; Arimidex/Tamoxifen","treatment type: ACx4 cycles",
                                                         "treatment type: ACx4; Tamoxifen","treatment type: ACx6; Tamoxifen","treatment type: ACx8 cycles; taxol",
                                                         "treatment type: CAF","treatment type: CAFx6 cycles","treatment type: ECx4; Tamoxifen",
                                                         "treatment type: Tam+AC","treatment type: Tam+ACx4 cycles","treatment type: Tam+ACx8;taxol",
                                                         "treatment type: Tam+CAFx2 cycles","treatment type: Tam+CAFx6 cycles"))
idNA = which(phenoData$characteristics_ch1.11   %in% c("treatment type: NA"))
idNone = which(phenoData$characteristics_ch1.11   %in% c("treatment type: none"))
covariate.df$anthra[idNA]=NA
covariate.df$anthra[idAntra]=T

idChemo = which(phenoData$characteristics_ch1.10 == "chemo? (0=no, 1=yes): 1")
idNA = which(phenoData$characteristics_ch1.10 == "chemo? (0=no, 1=yes): NA")
idNone = which(phenoData$characteristics_ch1.10 == "chemo? (0=no, 1=yes): 0")

covariate.df$chemo = NA
covariate.df$chemo[idChemo]=T
covariate.df$chemo[idNone]=F

idHm = which(phenoData$characteristics_ch1.9 == "characteristics: endocrine? (0=no, 1=yes): 1")
idNone = which(phenoData$characteristics_ch1.9 == "characteristics: endocrine? (0=no, 1=yes): 0")

covariate.df$ht = NA
covariate.df$ht[idHm]=T
covariate.df$ht[idNone]=F

covariate.df$NoTr =  covariate.df$chemo==F & covariate.df$ht==F

covariate.df$rfs.e = phenoData$characteristics_ch1.13=="dfs event (defined as any type of recurrence or death from breast cancer): 1"
covariate.df$rfs.t =  as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.12))*365

covariate.df$rdfs.e = phenoData$characteristics_ch1.15=="dmfs event (defined as distant metastasis or death from breast cancer): 1"
covariate.df$rdfs.t = as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.14))*365

covariate.df$rfs.e[(phenoData$characteristics_ch1.13=="dfs event (defined as any type of recurrence or death from breast cancer): NA" & phenoData$characteristics_ch1.15!="dmfs event (defined as distant metastasis or death from breast cancer): NA")]=covariate.df$rdfs.e[(phenoData$characteristics_ch1.13=="dfs event (defined as any type of recurrence or death from breast cancer): NA" & phenoData$characteristics_ch1.15!="dmfs event (defined as distant metastasis or death from breast cancer): NA")]
covariate.df$rfs.t[(phenoData$characteristics_ch1.13=="dfs event (defined as any type of recurrence or death from breast cancer): NA" & phenoData$characteristics_ch1.15!="dmfs event (defined as distant metastasis or death from breast cancer): NA")]=covariate.df$rdfs.t[(phenoData$characteristics_ch1.13=="dfs event (defined as any type of recurrence or death from breast cancer): NA" & phenoData$characteristics_ch1.15!="dmfs event (defined as distant metastasis or death from breast cancer): NA")]

covariate.df$os.e =  phenoData$characteristics_ch1.17=="dss event (defined as death from breast cancer): 1"
covariate.df$os.t = as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.16))*365
covariate.df$dss.e = phenoData$characteristics_ch1.17=="dss event (defined as death from breast cancer): 1"
covariate.df$pCR = NA

covariate.df$Herc = NA
covariate.df$Tam =covariate.df$ht
covariate.df$Tax= F
covariate.df$Tax[phenoData$characteristics_ch1.11   %in% c("treatment type: ACx8 cycles; taxol","treatment type: Tam+ACx8;taxol")]=T
covariate.df$Tax[phenoData$characteristics_ch1.11   %in% c("treatment type: NA")]=NA


load("../data/clinical/rmas/GSE45255.RData")
pData(eset)=covariate.df

featureData(eset)=featureData(rawsetL[[1]])
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/GSE45255/rmaData.RData")

