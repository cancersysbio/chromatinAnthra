library(GEOquery)

rawsetL = getGEO("GSE3494")
save(rawsetL,file = "../data/clinical/GSE3494_raw.RData")

load("../data/clinical/cells/GSE1456/Sweden_Clinical.RData")
load("../data/clinical/cells/GSE3494/uppsala-U133A.RData")

expOrig = exprs(rawsetL[[1]])
prExp = expOrig["208305_at",]
her2Exp = expOrig["216836_s_at",]
erExp = expOrig["205225_at",]


phenoData = read.table("../data/clinical/cells/GSE3494/clinic.txt",header = T,sep = "\t",stringsAsFactors = F)

## NO INFO TO SAY WHICH GET CMF or TAM or nothing, but none of them receive anthra (16846532)


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
thresholdPR = 4.294649
getThreshold(her2Exp)
thresholdher2 =9.601636
getThreshold(erExp)
thresholder =8.757161
pr_probe = prExp>thresholdPR
her2_probe = her2Exp>thresholdher2
er_probe = erExp>thresholder

probeTest = rownames(upp.x)[c(500:5500)]
exp_geo = exprs(rawsetL[[1]])[probeTest,]
exp_ext = upp.x[probeTest,]

dic_samples = NA
for(i in 1:ncol(exp_ext)){
  id_ext = colnames(upp.x)[i]
  vcor = apply(exp_geo,2,function(x) cor(x,exp_ext[,i]))
  idgeo = which.max(vcor)  
  maxgeo = max(vcor)
  
  dic_samples = rbind(dic_samples,c(id_ext,names(idgeo),maxgeo))
}
dic_samples = dic_samples[-c(1),]
SwedenClinical = SwedenClinical[ SwedenClinical$cohort=="Uppsala", ]

phenoData2 = cbind(phenoData,SwedenClinical[dic_samples[match(pData(rawsetL[[1]])$geo_accession,dic_samples[,2]),1],])

phenoData = phenoData2

covariate.df= data.frame(title=pData(rawsetL[[1]])$title, stringsAsFactors = F)
rownames(covariate.df)=colnames(rawsetL[[1]])
covariate.df$adjuvant = phenoData$treat.adjuv

covariate.df$neoadjuvant = F
covariate.df$adjNONE = !phenoData$treat.adjuv
covariate.df$adjNA = is.na(phenoData$treat.adjuv)


covariate.df$geoAc = pData(rawsetL[[1]])$geo_accession
covariate.df$cohort = "UPP"

covariate.df$age = phenoData$age.at.diagnosis
ER_1 = phenoData$ER.status=="ER+"
table(ER,ER_1)
table(er_probe,ER_1)

covariate.df$er = ER_1
PR_1 = phenoData$PgR.status=="PgR+"
table(PR,PR_1)
table(pr_probe,PR_1)

covariate.df$pr = PR


covariate.df$her2 = her2_probe
covariate.df$hist = NA
grade = as.numeric(gsub("[^0-9,.]", "",phenoData$Elston.histologic.grade)) 

covariate.df$grade = grade


size = as.numeric(phenoData$tumor.size..mm.)
lympNodePos = rep(NA,nrow(phenoData))
lympNodePos[which(phenoData$Lymph.node.status=="LN+")]=T
lympNodePos[which(phenoData$Lymph.node.status=="LN-")]=F



stage = NA

covariate.df$stage = stage
covariate.df$t_stage = ifelse(size>50,3,ifelse(size<20,1,2))
covariate.df$m_stage = NA
covariate.df$n_stage= NA
covariate.df$lympNodePos = lympNodePos
covariate.df$lymphNodeNum = NA
covariate.df$size= size
covariate.df$CX = NA
covariate.df$RX= NA
covariate.df$HX = NA

covariate.df$anthra = F

covariate.df$chemo = NA

covariate.df$ht = NA

covariate.df$NoTr =  covariate.df$adjuvant==F
covariate.df$NoTr[ is.na(covariate.df$adjuvant)]=NA


covariate.df$rfs.e = phenoData$DSS.EVENT..Disease.Specific.Survival.EVENT..1.death.from.breast.cancer..0.alive.or.censored..==1
covariate.df$rfs.e[is.na(phenoData$DSS.EVENT..Disease.Specific.Survival.EVENT..1.death.from.breast.cancer..0.alive.or.censored..)]=NA
covariate.df$rfs.t =  phenoData$DSS.TIME..Disease.Specific.Survival.Time.in.years.*365

covariate.df$rdfs.e =NA
covariate.df$rdfs.t =NA

covariate.df$os.e = covariate.df$rfs.e
covariate.df$os.t = covariate.df$rfs.t
covariate.df$dss.e = covariate.df$rfs.e
covariate.df$pCR = NA

covariate.df$Herc = NA
covariate.df$Tam =NA
covariate.df$Tax= NA

# process cells 

dir.create("../data/clinical/cells/GSE3494/GPL96/")
flist <- sapply(colnames(rawsetL[[1]]@assayData$exprs),function(x) list.files("../data/clinical/cells/GSE3494/",  paste(x,"*",sep = ""), full.names = TRUE))
file.copy(unlist(flist), "../data/clinical/cells/GSE3494/GPL96/")

dir.create("../data/clinical/cells/GSE3494/GPL97/")
flist <- sapply(colnames(rawsetL[[2]]@assayData$exprs),function(x) list.files("../data/clinical/cells/GSE3494/",  paste(x,"*",sep = ""), full.names = TRUE))
file.copy(unlist(flist), "../data/clinical/cells/GSE3494/GPL97/")


library(affy)

fns = list.celfiles("../data/clinical/cells/GSE3494/GPL96/",full.names=T)
celF = ReadAffy(filenames = fns)
eset <- rma(celF)
colnames(eset)=substr(colnames(eset),1,10)
pData(eset)=covariate.df

featureData(eset)=featureData(rawsetL[[1]])
eset = eset[which(rownames(eset) %in% eset@featureData@data$ID),]
dir.create("../data/clinical/rmas/GSE3494",recursive = T)
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/rmas/GSE3494/rmaData_96.RData")

fns = list.celfiles("../data/clinical/cells/GSE3494/GPL97/",full.names=T)
celF = ReadAffy(filenames = fns)
eset <- rma(celF)
colnames(eset)=substr(colnames(eset),1,8)
pData(eset)=cov

featureData(eset)=featureData(rawsetL[[2]])
#eset = eset[which(rownames(eset) %in% eset@featureData@data$ID),]
dir.create("../data/clinical/rmas/GSE3494",recursive = T)
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/GSE3494/rmaData_97.RData")
