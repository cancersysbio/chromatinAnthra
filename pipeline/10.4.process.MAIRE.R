library(GEOquery)

rawsetL = getGEO("GSE65194")
save(rawsetL,file = "./data/clinical//GSE65194_raw.RData")


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
thresholdPR = 2.261742
getThreshold(her2Exp)
thresholdher2 = 12.81415
getThreshold(erExp)
thresholder =4.476763
pr_probe = prExp>thresholdPR
her2_probe = her2Exp>thresholdher2
er_probe = erExp>thresholder

phenoData1 = read.table("/data2/data/BAF_cohort_files/cells/GSE65194/GSE65194_clinical_data_update.txt",header = T,sep = "\t",stringsAsFactors = F)[,1:40]
phenoData = GSE65194_UPC@phenoData@data
phenoData$ER_exp = ER
phenoData$PR_exp = PR
phenoData$Her2_exp = HER2
phenoData = base::merge(phenoData,phenoData1,by.x="geo_accession",by.y="GSE65194_ID")

#covariate.df = data.frame(title=rownames(GSE65194_UPC@phenoData@data), adjuvant =logical(), geoAc=character(), cohort=character(), age=numeric(),er=logical(),pr=logical(),her2=logical()
#                          ,hist=character(), grade=numeric(), stage=numeric(), t_stage=numeric()
#                          , m_stage=numeric(), n_stage=numeric(), lympNodePos=logical(), lymphNodeNum=numeric(), size=numeric(),
#                          CX=character(), RX=character(), HX=character(), anthra=logical(), rfs.e=logical(), rfs.t = numeric()
#                          , rdfs.e = logical(), rdfs.t = numeric(),os.t=numeric(),os.e=logical(),dss.e=logical(),pCR=logical(), stringsAsFactors = F)

covariate.df= data.frame(title=phenoData$title, stringsAsFactors = F)
rownames(covariate.df)=phenoData$geo_accession
covariate.df$adjuvant = T
covariate.df$neoadjuvant = F
covariate.df$adjNONE = is.na(phenoData$Chemotherapy)
covariate.df$adjNA = F

covariate.df$geoAc = phenoData$geo_accession
covariate.df$cohort = "MAIRE"

covariate.df$age = as.numeric(phenoData$Age.at.diagnosis.yrs.)
covariate.df$er = phenoData$ER_exp
covariate.df$pr = phenoData$PR_exp
covariate.df$her2 = phenoData$erbB2_POSITIVE
covariate.df$hist = phenoData$Histological.type.1

grade = NA
covariate.df$grade = NA

size = as.numeric(phenoData$Tumor.size.mm.)
t_stage = ifelse(size>50,3,ifelse(size<20,1,2))

ln = as.numeric(phenoData$Number.of.metastatic.lymph.nodes)
n_stage= NA
m_stage = NA

lymphNodeNum = as.numeric(phenoData$Number.of.metastatic.lymph.nodes)


stage = NA
covariate.df$stage = stage
covariate.df$t_stage = t_stage
covariate.df$m_stage = m_stage
covariate.df$n_stage= n_stage
covariate.df$lympNodePos = lymphNodeNum>0
covariate.df$lymphNodeNum = lymphNodeNum
covariate.df$size= size
covariate.df$CX = phenoData$CHEMO_DRUG_GROUP
covariate.df$RX= phenoData$Radiotherapy
covariate.df$HX = phenoData$HORMO_DRUG_GROUP
covariate.df$anthra = phenoData$CHEMO_DRUG_GROUP %in% c("Anthracyclines","Anthracyclines+Taxanes")
covariate.df$chemo = !is.na(phenoData$Chemotherapy)
covariate.df$ht = !is.na(phenoData$Hormone.therapy)
covariate.df$NoTr =  covariate.df$chemo==F

# RFS = local or distal or dead
covariate.df$rfs.e = phenoData$Local.recurrence=="yes" | phenoData$Metastasis=="yes" | phenoData$Death.event==1 
covariate.df$rfs.t =  pmin(as.numeric(phenoData$Local.recurrence.time.month.)*30.41,as.numeric(phenoData$Metastasis.time.month.)*30.41, as.numeric(phenoData$Death.event.time.months.)*30.41)
covariate.df$rdfs.e = phenoData$Metastasis=="yes" | phenoData$Death.event==1 
covariate.df$rdfs.t = pmin(as.numeric(phenoData$Metastasis.time.month.)*30.41, as.numeric(phenoData$Death.event.time.months.)*30.41)
                             
covariate.df$os.e = phenoData$Death.event==1 
covariate.df$os.t = as.numeric(phenoData$Death.event.time.months.)*30.41
covariate.df$dss.e = NA
covariate.df$pCR = NA

covariate.df$Herc = !is.na(phenoData$HERCEPTIN)
covariate.df$Tam =phenoData$HORMO_DRUG_GROUP %in% c("Aromatase inhibitors","Castration+Aromatase inhibitors","tamoxifen")
covariate.df$Tax= phenoData$CHEMO_DRUG_GROUP %in% c("Anthracyclines+Taxanes","Taxanes")



load("../data/clinical/rmas/GSE65194.RData")
pData(eset)=covariate.df

featureData(eset)=featureData(rawsetL[[1]])
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/GSE65194/rmaData.RData")