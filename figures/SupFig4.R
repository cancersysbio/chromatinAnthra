#Efig4

load("../data/complexes_v1.2_BRCA.RData")

load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")


gois2 = unique(mulvarrfs.res$gene[mulvarrfs.res$pval.expAnth<0.05])

probes_gois2 = mulvarrfs.res$probe[match(gois2,mulvarrfs.res$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

# forest plots
library(forestplot)
tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),signif(as.numeric(mulvarrfs.res$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res$probe)]),digits = 2) )


idhigh = which(unlist(contrast.res.l$lowExp[, "Contrast"])>unlist(contrast.res.l$highExp[, "Contrast"]))
idlow = which(unlist(contrast.res.l$lowExp[, "Contrast"])<unlist(contrast.res.l$highExp[, "Contrast"]))
idhigh = idhigh[order(unlist(contrast.res.l$highExp[idhigh, "Contrast"]))]
idlow = idlow[order(unlist(contrast.res.l$lowExp[idlow, "Contrast"]))]

pdf("EDFig4.pdf",width = 10,height = 20)
forestplot(tabletext[idhigh,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.2,
           mean = cbind(unlist(contrast.res.l$lowExp[idhigh, "Contrast"]),unlist(contrast.res.l$highExp[idhigh, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idhigh, "Lower"]), unlist(contrast.res.l$highExp[idhigh, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhigh, "Upper"]), unlist(contrast.res.l$highExp[idhigh, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlow,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.2,
           mean = cbind(unlist(contrast.res.l$lowExp[idlow, "Contrast"]),unlist(contrast.res.l$highExp[idlow, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idlow, "Lower"]), unlist(contrast.res.l$highExp[idlow, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlow, "Upper"]), unlist(contrast.res.l$highExp[idlow, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)

dev.off()

