library(GEOquery)

rawsetL = getGEO("GSE20685")
save(rawsetL,file = "./data/clinical//GSE20685_raw.RData")


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
thresholdPR = getThreshold(prExp)
thresholdPR = 6.033366
thresholdher2 = getThreshold(her2Exp)
thresholdher2 = 12.84642
thresholder = getThreshold(erExp)
thresholder =11.60443
pr_probe = prExp>thresholdPR
her2_probe = her2Exp>thresholdher2
er_probe = erExp>thresholder


phenoData = pData(rawsetL[[1]])

covariate.df= data.frame(title=phenoData$title, stringsAsFactors = F)
###############################################
rownames(covariate.df)=rownames(phenoData)
##############################################3
covariate.df$adjuvant = T
covariate.df$neoadjuvant = F
covariate.df$adjNONE = NA
covariate.df$adjNA = NA

covariate.df$geoAc = phenoData$geo_accession
covariate.df$cohort = "KAO"

covariate.df$age = as.numeric(gsub("[^0-9]", "", (phenoData$characteristics_ch1.2)))
covariate.df$er = er_probe
covariate.df$pr = pr_probe
covariate.df$her2 = her2_probe
covariate.df$hist = NA
covariate.df$grade = NA

id_1 = which(phenoData$characteristics_ch1.5=="event_metastasis: 1")
subtype = ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.7),as.character(phenoData$characteristics_ch1.6))
time_to_metastasis = ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.6),as.character(phenoData$characteristics_ch1.4))
t_stage = ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.8),as.character(phenoData$characteristics_ch1.7))
n_stage= ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.9),as.character(phenoData$characteristics_ch1.8))
m_stage = ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.10),as.character(phenoData$characteristics_ch1.9))
reg_rel = ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.11),as.character(phenoData$characteristics_ch1.10))

time_to_reg_rel = ifelse(reg_rel!="regional_relapse: 1",  as.character(phenoData$characteristics_ch1.4),
                         ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.12),as.character(phenoData$characteristics_ch1.11)))


adjuvant_chemotherapy = ifelse(reg_rel!="regional_relapse: 1",  ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.12),as.character(phenoData$characteristics_ch1.11)),
                               ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.13),as.character(phenoData$characteristics_ch1.12)))
cafVScmf =  ifelse(reg_rel!="regional_relapse: 1",  ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.13),as.character(phenoData$characteristics_ch1.12)),
                   ifelse(phenoData$characteristics_ch1.5=="event_metastasis: 1",as.character(phenoData$characteristics_ch1.14),as.character(phenoData$characteristics_ch1.13)))

covariate.df$adjuvant = adjuvant_chemotherapy=="adjuvant_chemotherapy: yes"
covariate.df$adjuvant[adjuvant_chemotherapy=="adjuvant_chemotherapy: unknown"] = NA
covariate.df$adjNONE = NA
covariate.df$adjNONE[adjuvant_chemotherapy=="adjuvant_chemotherapy: no"]=T
covariate.df$adjNONE[adjuvant_chemotherapy=="adjuvant_chemotherapy: yes"]=F
covariate.df$adjNA= F
covariate.df$adjNA[adjuvant_chemotherapy=="adjuvant_chemotherapy: unknown"]=T


lymphNodeNum = NA
size= NA


source("~/projects/BAF/tmn_stage.R")

stage = mapply(tmn_stage,t = as.numeric(gsub("[^0-9]", "", t_stage)),m = as.numeric(gsub("[^0-9]", "", m_stage)),n=as.numeric(gsub("[^0-9]", "", n_stage)))
covariate.df$stage = stage
covariate.df$t_stage = as.numeric(gsub("[^0-9]", "", t_stage))
covariate.df$m_stage = as.numeric(gsub("[^0-9]", "", m_stage))
covariate.df$n_stage= as.numeric(gsub("[^0-9]", "", n_stage))
covariate.df$lympNodePos = covariate.df$n_stage>0
covariate.df$lymphNodeNum = lymphNodeNum
covariate.df$size= size
covariate.df$CX = gsub("regimen (caf vs cmf): ","",cafVScmf,fixed = T)
covariate.df$RX= NA
covariate.df$HX = NA
covariate.df$anthra = NA
covariate.df$anthra[covariate.df$CX=="CAF"] = T
covariate.df$anthra[covariate.df$CX %in% c("CMF"," no")] = F
covariate.df$chemo = NA
covariate.df$chemo[covariate.df$CX %in% c("CAF","CMF","others")]=T
covariate.df$chemo[covariate.df$CX == "no"]=F
covariate.df$ht = NA
covariate.df$NoTr =  covariate.df$chemo==F


covariate.df$rfs.e = phenoData$characteristics_ch1.5=="event_metastasis: 1" | phenoData$characteristics_ch1.3=="event_death: 1"
covariate.df$rfs.t =  pmin(as.numeric(gsub("[^0-9,.]", "",time_to_metastasis))*365, as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.4))*365 )
covariate.df$rdfs.e = phenoData$characteristics_ch1.5=="event_metastasis: 1"
covariate.df$rdfs.t = as.numeric(gsub("[^0-9,.]", "",time_to_metastasis))*365
covariate.df$os.e = phenoData$characteristics_ch1.3=="event_death: 1"
covariate.df$os.t = as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.4))*365
covariate.df$dss.e = NA
covariate.df$pCR = NA
  
covariate.df$Herc = NA
covariate.df$Tam =NA
covariate.df$Tax= NA


load("../data/clinical/rmas/GSE20685.RData")
pData(eset)=covariate.df

featureData(eset)=featureData(rawsetL[[1]])
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/GSE20685/rmaData.RData")