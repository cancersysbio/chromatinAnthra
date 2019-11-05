load("../data/complexes_v1.2_BRCA.RData")

load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVStax_V2.RData")
mulvarrfs.res.aVt = mulvarrfs.res

load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSCMF_V2.RData")
mulvarrfs.res.aVc = mulvarrfs.res

load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")
mulvarrfs.res.aVn = mulvarrfs.res

load("../data/cell_line_dataset_VIPER_datasets_aac.RData")
gois2 = combine.pval$genes[which(combine.pval$combined<0.05)]
gois2= intersect(gois2,complexes$hgnc_symbol)

gois_avt = unique(mulvarrfs.res.aVt$gene[mulvarrfs.res.aVt$pval.expAnth<0.05])
gois_avc = unique(mulvarrfs.res.aVc$gene[mulvarrfs.res.aVc$pval.expAnth<0.05])
gois_avn = unique(mulvarrfs.res.aVn$gene[mulvarrfs.res.aVn$pval.expAnth<0.05])

load("../data/complexes_v1.2_BRCA.RData")
gois_avt = gois_avt[gois_avt %in% complexes$hgnc_symbol]
gois_avc = gois_avc[gois_avc %in% complexes$hgnc_symbol]
gois_avn = gois_avn[gois_avn %in% complexes$hgnc_symbol]



v = VennDiagram::venn.diagram(list("AnthraVSnoAnthra"=gois_avn,"AnthraVSTax"=gois_avt,"AnthraVSCMF"=gois_avc),
                              fill = viridis::viridis_pal()(3),filename=NULL)


pdf("EDFig7_e_venn.pdf",width = 4 ,height = 3)
grid::grid.draw(v)
dev.off()

viperNotInArray = setdiff(gois2,unique(mulvarrfs.res$gene))



genesInArray = mulvarrfs.res$gene
viperNotInArray = setdiff(gois2,genesInArray)

library(viridis)



# non anthra vs anthra
v = VennDiagram::venn.diagram(list("Viper"=gois2,"Metacohort"=gois_avn,"not in array"=viperNotInArray),
                              fill = viridis_pal()(3),filename=NULL)
pdf("EDFig7_a_venn_anthra_nonAntra.pdf",width = 4,height = 3)
grid::grid.draw(v)
dev.off()

# cmf vs anthra
pdf("EDFig7_b_venn_anthra_cmf.pdf")
v2 = VennDiagram::venn.diagram(list("viper"=gois2,"clinic"=gois_avc,"notInArray"=viperNotInArray),
                               fill =viridis_pal()(3),filename=NULL)
grid::grid.draw(v2)
dev.off()



# tax vs anthra
v = VennDiagram::venn.diagram(list("Viper"=gois2,"Metacohort"=gois_avt,"not in array"=viperNotInArray),
                              fill = viridis_pal()(3),filename=NULL)
pdf("EDFig7_c_venn_anthra_tax.pdf",width = 4,height = 3)
grid::grid.draw(v)
dev.off()


# plot forest
library(forestplot)


# anthra vs non anthra
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")


probes_gois2 = mulvarrfs.res$probe[match(gois2,mulvarrfs.res$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),
                  signif(as.numeric(mulvarrfs.res.aVn$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res.aVn$probe)]),digits = 2)
)
#View(mulvarrfs.res.aVn[which(as.numeric(as.character(mulvarrfs.res.aVn$pval.expAnth))<0.05&as.numeric(as.character(mulvarrfs.res.aVn$odds.expAnth))>1&mulvarrfs.res.aVn$gene %in% tabletext[,1]),])

idhighExp = c(19,20,3,23,29,10) #c(12,13,2,16)
idlowexp =   c(6,26,17,7,28,11) #c(17,18,9,7,19)
pdf("Fig3c.pdf",width = 10,height = 10)
forestplot(tabletext[idhighExp,-c(2)],
           legend = c("lowExp", "highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Contrast"]), unlist(contrast.res.l$highExp[idhighExp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Lower"]),  unlist(contrast.res.l$highExp[idhighExp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Upper"]),  unlist(contrast.res.l$highExp[idhighExp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlowexp,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Contrast"]), unlist(contrast.res.l$highExp[idlowexp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Lower"]),  unlist(contrast.res.l$highExp[idlowexp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Upper"]),  unlist(contrast.res.l$highExp[idlowexp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)
dev.off()

# anthra vs cmf
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSCMF_V2.RData")


probes_gois2 = mulvarrfs.res.aVc$probe[match(gois2,mulvarrfs.res.aVc$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),
                  signif(as.numeric(mulvarrfs.res.aVc$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res.aVc$probe)]),digits = 2) )
#View(mulvarrfs.res.aVc[which(as.numeric(as.character(mulvarrfs.res.aVc$pval.expAnth))<0.05&as.numeric(as.character(mulvarrfs.res.aVc$odds.expAnth))<1&mulvarrfs.res.aVc$gene %in% tabletext[,1]),])

idhighExp = c(17,22,3,20,8) 
idlowexp =  c(14,25,5,19,27,9) 
pdf("Fig3d.pdf",width = 10,height = 10)
forestplot(tabletext[idhighExp,-c(2)],
           legend = c("lowExp", "highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Contrast"]), unlist(contrast.res.l$highExp[idhighExp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Lower"]),  unlist(contrast.res.l$highExp[idhighExp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Upper"]),  unlist(contrast.res.l$highExp[idhighExp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlowexp,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Contrast"]), unlist(contrast.res.l$highExp[idlowexp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Lower"]),  unlist(contrast.res.l$highExp[idlowexp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Upper"]),  unlist(contrast.res.l$highExp[idlowexp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)
dev.off()

# anthra vs tax
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVStax_V2.RData")


probes_gois2 = mulvarrfs.res.aVt$probe[match(gois2,mulvarrfs.res.aVt$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),
                  signif(as.numeric(mulvarrfs.res.aVt$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res.aVt$probe)]),digits = 2) )


idhighExp = c(17,23) 
idlowexp =  c(11,31)
pdf("Fig3e.pdf",width = 10,height = 10)
forestplot(tabletext[idhighExp,-c(2)],
           legend = c("lowExp", "highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Contrast"]), unlist(contrast.res.l$highExp[idhighExp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Lower"]),  unlist(contrast.res.l$highExp[idhighExp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhighExp, "Upper"]),  unlist(contrast.res.l$highExp[idhighExp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlowexp,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Contrast"]), unlist(contrast.res.l$highExp[idlowexp, "Contrast"])),
           lower = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Lower"]),  unlist(contrast.res.l$highExp[idlowexp, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlowexp, "Upper"]),  unlist(contrast.res.l$highExp[idlowexp, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)
dev.off()


# Her2+
load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_MKI67_anthraVSnoAnthra_V2.RData")
contrast.res =contrast.res.her2P

probes_gois2 = mulvarrfs.res.her2P$probe[match(gois2,mulvarrfs.res.her2P$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),signif(as.numeric(mulvarrfs.res.her2P$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res.her2P$probe)]),digits = 2) )

idhigh = c(24,20,3,14,29,11) 
idlow =   c(8,27,18,23,19,12)

pdf("EDFig8_her2P.pdf",width = 10,height = 10)
forestplot(tabletext[idhigh,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhigh, "Contrast"]),unlist(contrast.res.l$highExp[idhigh, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idhigh, "Lower"]), unlist(contrast.res.l$highExp[idhigh, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhigh, "Upper"]), unlist(contrast.res.l$highExp[idhigh, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlow,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlow, "Contrast"]),unlist(contrast.res.l$highExp[idlow, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idlow, "Lower"]), unlist(contrast.res.l$highExp[idlow, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlow, "Upper"]), unlist(contrast.res.l$highExp[idlow, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)

dev.off()


# ER+Her2-
load("../data/metacohort_assoc_results_age_cohort_er_pr_her2_tsaeg_ln_anthraVSnoAnthra_V2_ERP_withTam.RData")
load("../data/metacohort_assoc_contrasts_age_cohort_er_pr_her2_tsaeg_ln_anthraVSnoAnthra_V2_ERP_withTam.RData")
contrast.res =contrast.res.erP

probes_gois2 = mulvarrfs.res.erP$probe[match(gois2,mulvarrfs.res.erP$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),697)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),signif(as.numeric(mulvarrfs.res.erP$pval.expAnth[match(as.character(contrast.res.l$highExp$probe),mulvarrfs.res.erP$probe)]),digits = 2) )

idhigh = c(20,24,3,22,29,11) 
idlow =   c(8,27,17,9,18,12)

pdf("EDFig8_ERPHer2N.pdf",width = 10,height = 10)
forestplot(tabletext[idhigh,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhigh, "Contrast"]),unlist(contrast.res.l$highExp[idhigh, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idhigh, "Lower"]), unlist(contrast.res.l$highExp[idhigh, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhigh, "Upper"]), unlist(contrast.res.l$highExp[idhigh, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlow,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlow, "Contrast"]),unlist(contrast.res.l$highExp[idlow, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idlow, "Lower"]), unlist(contrast.res.l$highExp[idlow, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlow, "Upper"]), unlist(contrast.res.l$highExp[idlow, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)

dev.off()


# TNBC
contrast.res =contrast.res.TNBC

probes_gois2 =  mulvarrfs.res.TNBC$probe[match(gois2, mulvarrfs.res.TNBC$gene)]

contrast.res$desc = rep.int(c("lowExp","highExp"),552)
g = as.factor(contrast.res$desc)
contrast.res.l = split(contrast.res[contrast.res$probe %in% probes_gois2,],g[contrast.res$probe %in% probes_gois2])

tabletext = cbind(as.character(contrast.res.l$highExp$gene),as.character(contrast.res.l$highExp$probe),signif(as.numeric( mulvarrfs.res.TNBC$pval.expAnth[match(as.character(contrast.res.l$highExp$probe), mulvarrfs.res.TNBC$probe)]),digits = 2) )

idhigh = c(18,24,3,21,29,9) 
idlow =   c(5,26,14,6,28,10)

pdf("EDFig8_TNBC.pdf",width = 10,height = 10)
forestplot(tabletext[idhigh,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idhigh, "Contrast"]),unlist(contrast.res.l$highExp[idhigh, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idhigh, "Lower"]), unlist(contrast.res.l$highExp[idhigh, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idhigh, "Upper"]), unlist(contrast.res.l$highExp[idhigh, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = F
)

forestplot(tabletext[idlow,-c(2)],
           legend = c("lowExp","highExp"),
           fn.ci_norm =fpDrawNormalCI,boxsize = 0.1,
           mean = cbind(unlist(contrast.res.l$lowExp[idlow, "Contrast"]),unlist(contrast.res.l$highExp[idlow, "Contrast"])), 
           lower = cbind(unlist(contrast.res.l$lowExp[idlow, "Lower"]), unlist(contrast.res.l$highExp[idlow, "Lower"])),
           upper = cbind(unlist(contrast.res.l$lowExp[idlow, "Upper"]), unlist(contrast.res.l$highExp[idlow, "Upper"])),
           col=fpColors(box=viridis_pal()(3)[-c(3)]),
           xlab="OS hazard",new_page = T
)

dev.off()