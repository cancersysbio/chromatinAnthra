#figure 3

#a
load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_anthraVSnoAnthra_V2_tritorax.RData")
complex_trit = read.table("../data/trithoraxVSPRC.tsv.csv",header = T,sep="\t",stringsAsFactors = F)

mulvarrfs.res$group = complex_trit$Group[match(mulvarrfs.res$entrez,complex_trit$Entrez)]
mulvarrfs.res$complex = complex_trit$Complex[match(mulvarrfs.res$entrez,complex_trit$Entrez)]
mulvarrfs.res$complex2 = NA
mulvarrfs.res$complex2[mulvarrfs.res$complex %in% c("CanonicalPRC1","CorePRC1","NonCanonicalPRC1")]="PRC1"
mulvarrfs.res$complex2[mulvarrfs.res$complex %in% c("CorePRC2","PRC2AccesoryProt")]="PRC2"
mulvarrfs.res$complex2[mulvarrfs.res$complex %in% c("CorePR-DUB","PR-DUBAccesory")]="PR-DUB"
mulvarrfs.res$complex2[mulvarrfs.res$complex %in% c("SWI/SNF")]="SWI/SNF"
mulvarrfs.res$complex2[mulvarrfs.res$complex %in% c("ASH1","CoreCompass","MLL1/Compass",
                                                    "MLL3/Compass","SET1/Compass","")]="COMPASS"


mulvarrfs.res_trit_pcr = mulvarrfs.res[which(mulvarrfs.res$entrez %in% na.omit(complex_trit$Entrez[which(complex_trit$Group=="PcG") ] )),]
mulvarrfs.res_trit_pcr = mulvarrfs.res_trit_pcr[which(!duplicated(mulvarrfs.res_trit_pcr$entrez)),] 
mulvarrfs.res_trit_trit = mulvarrfs.res[which(mulvarrfs.res$entrez %in% na.omit(complex_trit$Entrez[which(complex_trit$Group=="TrxG") ] )),]
mulvarrfs.res_trit_trit=mulvarrfs.res_trit_trit[which(!duplicated(mulvarrfs.res_trit_trit$entrez)),]


library(dplyr)
complex.est.pcr= ddply(mulvarrfs.res_trit_pcr,.(complex2),summarise,pval = survcomp::combine.test(as.numeric(pval.expAnth)))

complex.est.trit = ddply(mulvarrfs.res_trit_trit,.(complex2),summarise,pval = survcomp::combine.test(as.numeric(pval.expAnth)))


complex.est.pcr = cbind(complex.est.pcr,ddply(mulvarrfs.res_trit_pcr,.(complex2),summarise,
                                              coef = survcomp::combine.est(as.numeric(coef.expAntha),as.numeric(error.odds.expAnth))[1],
                                              error = survcomp::combine.est(as.numeric(coef.expAntha),as.numeric(error.odds.expAnth))[2])[,2:3])

complex.est.trit = cbind(complex.est.trit,ddply(mulvarrfs.res_trit_trit,.(complex2),summarise,
                                                coef = survcomp::combine.est(as.numeric(coef.expAntha),as.numeric(error.odds.expAnth))[1],
                                                error = survcomp::combine.est(as.numeric(coef.expAntha),as.numeric(error.odds.expAnth))[2])[,2:3])



complex.est.pcr$low = apply(complex.est.pcr ,1,function(x) (x$coef)+qnorm(0.025, lower.tail=TRUE,log.p = F)*x$error)
complex.est.pcr$high = apply(complex.est.pcr ,1,function(x) (x$coef)+qnorm(0.025, lower.tail=F)*x$error)
complex.est.trit$low = apply(complex.est.trit ,1,function(x) (x$coef)+qnorm(0.025, lower.tail=TRUE)*x$error)
complex.est.trit$high = apply(complex.est.trit ,1,function(x) (x$coef)+qnorm(0.025, lower.tail=F)*x$error)


complex.est.pcr$pval2 = signif(complex.est.pcr$pval,2)
complex.est.trit$pval2 = signif(complex.est.trit$pval,2)

complex.est.pcr$ci = paste("(",signif(complex.est.pcr$low,3)," ",signif(complex.est.pcr$high,3),")",sep="")
complex.est.trit$ci = paste("(",signif(complex.est.trit$low,3)," ",signif(complex.est.trit$high,3),")",sep="")



pdf("Fig3a.pdf")

forestplot.surv(labeltext = complex.est.pcr[,c("complex2","ci"),drop=F],mean=(as.numeric(complex.est.pcr$coef)),
                lower=as.numeric(complex.est.pcr$low),upper = as.numeric(complex.est.pcr$high),xlog = F,zero=0)
forestplot.surv(labeltext = complex.est.trit[,c("complex2","ci"),drop=F],mean=(as.numeric(complex.est.trit$coef)),
                lower=as.numeric(complex.est.trit$low),upper = as.numeric(complex.est.trit$high),xlog=F,zero=0)


dev.off()


# figure ED6

mulvarrfs.res_trit_trit$pval2 = signif(as.numeric(mulvarrfs.res_trit_trit$pval.expAnth),2)
mulvarrfs.res_trit_pcr$pval2 = signif(as.numeric(mulvarrfs.res_trit_pcr$pval.expAnth),2)

mulvarrfs.res_trit_trit$ci =    paste("(",signif(as.numeric(mulvarrfs.res_trit_trit$conf.l.expAntra),3)," ",signif(as.numeric(mulvarrfs.res_trit_trit$conf.h.expAntra),3),")",sep="")
mulvarrfs.res_trit_pcr$ci =    paste("(",signif(as.numeric(mulvarrfs.res_trit_pcr$conf.l.expAntra),3)," ",signif(as.numeric(mulvarrfs.res_trit_pcr$conf.h.expAntra),3),")",sep="")

pdf("EDFig6.pdf")

# plot forest for each complex
# SWI/SNF

comb.stat = survcomp::combine.est(as.numeric(mulvarrfs.res_trit_trit$coef.expAntha[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF")]),
                                  as.numeric(mulvarrfs.res_trit_trit$error.odds.expAnth[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF")]))


comb.stat$ci.l = comb.stat$estimate+qnorm(0.025, lower.tail=T)*comb.stat$se
comb.stat$ci.h = comb.stat$estimate+qnorm(0.025, lower.tail=F)*comb.stat$se

forestplot.surv(labeltext = rbind(mulvarrfs.res_trit_trit[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF"),c("gene","pval2","ci"),drop=F],c("Combined","",paste("(",signif(comb.stat$ci.l,3)," ",signif(comb.stat$ci.h,3),")",sep="")) )
                ,mean=c(as.numeric(mulvarrfs.res_trit_trit$coef.expAntha[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF")]),comb.stat$estimate,comb.stat$estimate),
                lower=c(as.numeric(mulvarrfs.res_trit_trit$conf.l.expAntra[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF")]),comb.stat$ci.l,comb.stat$ci.l),
                upper = c(as.numeric(mulvarrfs.res_trit_trit$conf.h.expAntra[which(mulvarrfs.res_trit_trit$complex2=="SWI/SNF")]),comb.stat$ci.h,comb.stat$ci.h),
                is.summary = c(rep(F,23),T),xlog=F,zero=0)


#compass

comb.stat = survcomp::combine.est(as.numeric(mulvarrfs.res_trit_trit$coef.expAntha[which(mulvarrfs.res_trit_trit$complex2=="COMPASS")]),
                                  as.numeric(mulvarrfs.res_trit_trit$error.odds.expAnth[which(mulvarrfs.res_trit_trit$complex2=="COMPASS")]))

comb.stat$ci.l = comb.stat$estimate+qnorm(0.025, lower.tail=T)*comb.stat$se
comb.stat$ci.h = comb.stat$estimate+qnorm(0.025, lower.tail=F)*comb.stat$se

forestplot.surv(labeltext = rbind(mulvarrfs.res_trit_trit[which(mulvarrfs.res_trit_trit$complex2=="COMPASS"),c("gene","pval2","ci"),drop=F],c("Combined","",paste("(",signif(comb.stat$ci.l,3)," ",signif(comb.stat$ci.h,3),")",sep="") ) )
                ,mean=c(as.numeric(mulvarrfs.res_trit_trit$coef.expAntha[which(mulvarrfs.res_trit_trit$complex2=="COMPASS")]),comb.stat$estimate,comb.stat$estimate),
                lower=c(as.numeric(mulvarrfs.res_trit_trit$conf.l.expAntra[which(mulvarrfs.res_trit_trit$complex2=="COMPASS")]),comb.stat$ci.l,comb.stat$ci.l),
                upper = c(as.numeric(mulvarrfs.res_trit_trit$conf.h.expAntra[which(mulvarrfs.res_trit_trit$complex2=="COMPASS")]),comb.stat$ci.h,comb.stat$ci.h),
                is.summary = c(rep(F,17),T),xlog=F,zero=0)

# prc1

comb.stat = survcomp::combine.est(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PRC1")]),
                                  as.numeric(mulvarrfs.res_trit_pcr$error.odds.expAnth[which(mulvarrfs.res_trit_pcr$complex2=="PRC1")]))

comb.stat$ci.l = comb.stat$estimate+qnorm(0.025, lower.tail=T)*comb.stat$se
comb.stat$ci.h = comb.stat$estimate+qnorm(0.025, lower.tail=F)*comb.stat$se

forestplot.surv(labeltext = rbind(mulvarrfs.res_trit_pcr[which(mulvarrfs.res_trit_pcr$complex2=="PRC1"),c("gene","pval2","ci"),drop=F],c("Combined","",paste("(",signif(comb.stat$ci.l,3)," ",signif(comb.stat$ci.h,3),")",sep="") ) )
                ,mean=c(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PRC1")]),comb.stat$estimate,comb.stat$estimate),
                lower=c(as.numeric(mulvarrfs.res_trit_pcr$conf.l.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PRC1")]),comb.stat$ci.l,comb.stat$ci.l),
                upper = c(as.numeric(mulvarrfs.res_trit_pcr$conf.h.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PRC1")]),comb.stat$ci.h,comb.stat$ci.h),
                is.summary = c(rep(F,17),T),xlog=F,zero=0)



# prc2

comb.stat = survcomp::combine.est(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PRC2")]),
                                  as.numeric(mulvarrfs.res_trit_pcr$error.odds.expAnth[which(mulvarrfs.res_trit_pcr$complex2=="PRC2")]))

comb.stat$ci.l = comb.stat$estimate+qnorm(0.025, lower.tail=T)*comb.stat$se
comb.stat$ci.h = comb.stat$estimate+qnorm(0.025, lower.tail=F)*comb.stat$se

forestplot.surv(labeltext = rbind(mulvarrfs.res_trit_pcr[which(mulvarrfs.res_trit_pcr$complex2=="PRC2"),c("gene","pval2","ci"),drop=F],c("Combined","",paste("(",signif(comb.stat$ci.l,3)," ",signif(comb.stat$ci.h,3),")",sep="") ) )
                ,mean=c(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PRC2")]),comb.stat$estimate,comb.stat$estimate),
                lower=c(as.numeric(mulvarrfs.res_trit_pcr$conf.l.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PRC2")]),comb.stat$ci.l,comb.stat$ci.l),
                upper = c(as.numeric(mulvarrfs.res_trit_pcr$conf.h.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PRC2")]),comb.stat$ci.h,comb.stat$ci.h),
                is.summary = c(rep(F,9),T),xlog=F,zero=0)

# PR-DUB

comb.stat = survcomp::combine.est(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB")]),
                                  as.numeric(mulvarrfs.res_trit_pcr$error.odds.expAnth[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB")]))

comb.stat$ci.l = comb.stat$estimate+qnorm(0.025, lower.tail=T)*comb.stat$se
comb.stat$ci.h = comb.stat$estimate+qnorm(0.025, lower.tail=F)*comb.stat$se

forestplot.surv(labeltext = rbind(mulvarrfs.res_trit_pcr[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB"),c("gene","pval2","ci"),drop=F],c("Combined","",paste("(",signif(comb.stat$ci.l,3)," ",signif(comb.stat$ci.h,3),")",sep="") ))
                ,mean=c(as.numeric(mulvarrfs.res_trit_pcr$coef.expAntha[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB")]),comb.stat$estimate,comb.stat$estimate),
                lower=c(as.numeric(mulvarrfs.res_trit_pcr$conf.l.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB")]),comb.stat$ci.l,comb.stat$ci.l),
                upper = c(as.numeric(mulvarrfs.res_trit_pcr$conf.h.expAntra[which(mulvarrfs.res_trit_pcr$complex2=="PR-DUB")]),comb.stat$ci.h,comb.stat$ci.h),
                is.summary = c(rep(F,6),T),xlog=F,zero=0)

dev.off()

