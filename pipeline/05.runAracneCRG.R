library(org.Hs.eg.db)
library(RTN)

load("../data/complexes_v1.2_BRCA.RData")



load("../data/vst.RData")


des_raw = select(org.Hs.eg.db, rownames(exp.T.vst), c("ENTREZID","GENENAME","ENSEMBL"), "ENSEMBL")
exp.T.vst.e = exp.T.vst
exp.N.vst.e = exp.N.vst
Des$entrez = des_raw$ENTREZID[match(Des$ensembl_gene_id,des_raw$ENSEMBL)]
rownames(exp.T.vst.e)=Des$external_gene_name
rownames(exp.N.vst.e)=Des$external_gene_name
idna = which(is.na(Des$entrez))
#idna = which(rownames(exp.T.vst.e)=="")
Des.e = Des
exp.T.vst.e=exp.T.vst.e[-idna,]
exp.N.vst.e=exp.N.vst.e[-idna,]
Des.e = Des.e[-idna,]
iddup = duplicated(rownames(exp.T.vst.e))
exp.T.vst.e= exp.T.vst.e[which(!iddup),]
exp.N.vst.e= exp.N.vst.e[which(!iddup),]
Des.e = Des.e[which(!iddup),]
geneNames = unique(Des.e$external_gene_name[Des.e$entrez %in% unique(complexes$entrezgene)])
names(geneNames)=geneNames
tumor.TNI <- new("TNI", gexp=exp.T.vst.e, regulatoryElements=geneNames)
tumor.rtni <- tni.preprocess(tumor.TNI)


tumor.rtni<-tni.permutation(tumor.rtni)
tumor.rtni<-tni.bootstrap(tumor.rtni)
tumor.rtni<-tni.dpi.filter(tumor.rtni)

save(tumor.rtni,file="../data/allTumor_rtni.RData")

#aracne
load("../data/allTumor_rtni.RData")


tpc2regulon = function(net,filterTFs=NULL,pvalueCutoff= 0.05){
    #net = normal.rtni
    
    if (is.null(filterTFs)) {TFs <- colnames(net@results$tn.dpi)}
    if (!is.null(filterTFs)) {TFs <- filterTFs; TFs <- intersect(TFs, colnames(net@results$tn.dpi))}
    
    
    mode1 = sign(net@results$tn.dpi[,TFs])
    
    scores <- net@results$tn.dpi[,TFs]
    qvalues = net@results$adjpv[,TFs]
    
    
    aracne <- list()
    for (tf in TFs) {
      reg <- qvalues[,tf]
      which.reg <- reg <= pvalueCutoff
      which.mode = mode1[,tf]!=0
      likelihood <- scores[which.reg&which.mode,match(tf, colnames(scores))]
      tfmode <- mode1[which.reg&which.mode,match(tf, colnames(mode1))]
      aracne[[tf]] <- list("tfmode"=tfmode, "likelihood"=likelihood)
    }
    
    # removing missing data from the aracne regulon
    aracne <- aracne[names(aracne) != "NA"]
    aracne <- lapply(aracne, function(x) {
      filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
      x$tfmode <- x$tfmode[filtro]
      x$likelihood <- x$likelihood[filtro]
      return(x)
    })
    
    regul <- aracne[sapply(aracne, function(x) length(names(x$tfmode)))>0]
    class(regul) <- "regulon"
    
    return(regul)
  }


regul2.T = tpc2regulon(tumor.rtni)
chromatin_regulon = regul2.T
save(chromatin_regulon,file = "../data/chromatin_regulon.RData")





#Supfig2



chromatin_degrees = unlist(lapply(chromatin_regulon,function(x) length(x$tfmode)))
# get histogram of regulon
toPlot = data.frame(gene=names(chromatin_degrees),degree=chromatin_degrees,stringsAsFactors = F)

ggpubr::gghistogram(toPlot,x="degree",y="..count..",xlab = "Degree",rug = TRUE ,add="median",color=viridis::viridis_pal()(3)[2],fill=viridis::viridis_pal()(3)[2] )
ggsave("../plots/SupFig2_histogram_CRG_regulon.pdf",width = 7,height = 4)
# get barplot regulon

ggpubr::ggbarplot(toPlot[toPlot$degree>500,],x="gene",y="degree",sort.val = "asc",
                  color=viridis::viridis_pal()(3)[1],fill=viridis::viridis_pal()(3)[1])+theme(axis.text.x = element_text(angle = 35,vjust=0.65))
ggsave("../plots/SupFig2_barplot_CRG_regulon.pdf",width = 7,height = 4)
