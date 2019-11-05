library(rms)
library(viridis)
load("../data/clinical_metacohort_anthra_vs_nonAnthra.RData")
load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")


rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}


plot_geneInfo = function(mergedCOMBAT,mulvarrfs.res,pheno.df,genename="KDM4B",label){
  
  #KDM4B
  pdf(paste0("surv_",label,"_",genename,".pdf"),width = 8,height = 7)
  #i=which(mergedCOMBAT@featureData@data$`Gene Symbol`==genename)
  #probes_i = mergedCOMBAT@featureData@data$ID[i]
  probes= mulvarrfs.res$probe[which(mulvarrfs.res$gene==genename)]
  
  X.b = cbind(pheno.df[,setdiff(colnames(pheno.df),"tmain")],exp_c=rescale(exprs(mergedCOMBAT)[probes[1],rownames(pheno.df) ]))
  d <- datadist(X.b)
  #options(datadist='dd')
  assign('dd', datadist(X.b))
  options(datadist="dd")
  cxmodel2 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+er+pr+her2+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b,x = T,y=T,surv = T)
  
  print(cxmodel2)
  print(anova(cxmodel2))
  
  min_lim = median(X.b$exp_c)-sd(X.b$exp_c)
  max_lim = median(X.b$exp_c)+sd(X.b$exp_c)
  
  # low antrha green / low non anthra violet
  # high anthra yellow / high non anthra blue
  survplot(cxmodel2, anthra,exp_c=min_lim,er=F,pr=F,her2=F,conf.int = T,
           adj.subtitle=F,col=viridis_pal()(4)[c(1,3)],xlim=c(0,13),lwd=3,lty=1,label.curves=F)
  survplot(cxmodel2, anthra,exp_c=max_lim,er=F,pr=FALSE,her2=F,conf.int = T,
           adj.subtitle=F,col=viridis_pal()(4)[c(2,4)],xlim=c(0,13),lwd=3,lty=1,label.curves=F)
  
  X.b$exp_chemo =NA
  new_exp = survminer::surv_categorize( survminer::surv_cutpoint(X.b,time="rfs.t",event="rfs.e",variables = "exp_c"))
  plot(survminer::surv_cutpoint(X.b,time="rfs.t",event="rfs.e",variables = "exp_c"))
  X.b$exp_chemo[new_exp$exp_c=="low" & X.b$anthra=="TRUE"]="Low expression and anthracycline"
  X.b$exp_chemo[new_exp$exp_c=="low" & X.b$anthra=="FALSE"]="Low expression and non anthracycline"
  X.b$exp_chemo[new_exp$exp_c=="high" & X.b$anthra=="TRUE"]="High expression and anthracycline"
  X.b$exp_chemo[new_exp$exp_c=="high" & X.b$anthra=="FALSE"]="High expression and non anthracycline"
  
  survcomp::km.coxph.plot(Surv(rfs.t,rfs.e)~exp_chemo,data.s = X.b,
                          x.label="time (years)",y.label="OS",main.title="",show.n.risk = T,
                          n.risk.step=3,leg.text = levels(factor(X.b$exp_chemo)),
                          .col= viridis_pal()(4)[c(4,2,3,1)],xlim=c(0,13),.lwd=3)
  
  
  
  
  
  X.b$exp_chemo_2 = X.b$exp_chemo
  X.b$exp_chemo_2[X.b$exp_chemo %in% c("High expression and anthracycline","High expression and non anthracycline")]=NA
  survcomp::km.coxph.plot(Surv(rfs.t,rfs.e)~exp_chemo_2,data.s = X.b,
                          x.label="time (years)",y.label="OS",main.title=" ",show.n.risk = T,
                          n.risk.step=3,leg.text = levels(factor(X.b$exp_chemo_2)),
                          .col= viridis_pal()(4)[c(3,1)],xlim=c(0,13),.lwd=3)
  
  
  
  X.b$exp_chemo_3 = X.b$exp_chemo
  X.b$exp_chemo_3[!X.b$exp_chemo %in% c("High expression and anthracycline","High expression and non anthracycline")]=NA
  survcomp::km.coxph.plot(Surv(rfs.t,rfs.e)~exp_chemo_3,data.s = X.b,
                          x.label="time (years)",y.label="OS",main.title=" ",show.n.risk = T,
                          n.risk.step=3,leg.text = levels(factor(X.b$exp_chemo_3)),
                          .col= viridis_pal()(4)[c(4,2)],xlim=c(0,13),.lwd=3)
  #table(X.b$exp_chemo_3,X.b$cohort)
  #table(X.b$exp_chemo_3,X.b$X20)
  
  cxmodel3 = cph(Surv(rfs.t,rfs.e)~exp_c*anthra+rcs(age,3)+er+pr+her2+lymphNodePos+t.stage+cohort+MKI67_4,data = X.b,x = T,y=T,surv = T)
  hazards3 = predict(cxmodel2,X.b,se.fit = T,center.terms=T)

  
  p = ggplot(Predict(cxmodel3,exp_c,anthra,er=T,pr=F,her2=F,t.stage=2,lymphNodePos=F,cohort="KAO"))+theme_bw()+scale_color_viridis(discrete = T)+
    geom_point(data =X.b,aes(x=exp_c,y=hazards3$linear.predictors))+
    xlim(dd$limits["Low:prediction","exp_c"],dd$limits["High:prediction","exp_c"])
  
  print(p)
  
  dev.off()
  
  
}

plot_geneInfo(clinical_metacohort,mulvarrfs.res,pheno.df,"KDM4B","anthra_vs_nonAnthra_v2")
plot_geneInfo(clinical_metacohort,mulvarrfs.res,pheno.df,"KAT6B","anthra_vs_nonAnthra_v2")
plot_geneInfo(clinical_metacohort,mulvarrfs.res,pheno.df,"BCL11A","anthra_vs_nonAnthra_v2")


