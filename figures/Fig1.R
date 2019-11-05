# Figure 1 c
load("../data/cellLineDoxoSignature.RData")
load("../data/chromatin_regulon.RData")

# to fig1b,c
doxoSig.n = doxoSig
nullmodel.n =nullmodel
chromatin_regulon.n = chromatin_regulon
library(viper)
mrs.T <- msviper(doxoSig.n[["statistic"]], chromatin_regulon.n, nullmodel.n, verbose =T,pleiotropy = F)
p1 <- sort(mrs.T$es$p.value)
MRs <- names(p1)[p1 < 0.1]
plot(mrs.T, length(MRs))


