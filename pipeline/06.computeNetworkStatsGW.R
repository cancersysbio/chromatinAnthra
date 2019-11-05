
library(igraph)


aracneFile = "../data/ARACNE_all_res/network.txt"

metrics_weigthed_ARACNE2 = function(aracneFile = "network.txt",folder=".",label="",genelist=NULL){
  #http://bioinformatics.oxfordjournals.org/content/26/2/205.long
  
  
  
  
  tumor = read.table(aracneFile,header = T,sep = "\t",stringsAsFactors = F)
  tumor$MI[tumor$ MI<=0]=1E-5
  if(!is.null(genelist)){
    tumor = tumor[which(tumor$Regulator %in% genelist & tumor$Target %in% genelist),]
  }
  gtumor = graph.data.frame(tumor)
  E(gtumor)$weight = E(gtumor)$MI
  
  
  geneList = unique(c(tumor$Regulator,tumor$Target))  
  
  
  graphStats = data.frame(genes=geneList,stringsAsFactors = F)
  rownames(graphStats) = graphStats$genes
  print("degree")
  stat1 = graph.strength(gtumor,mode="all",loops=F)
  graphStats$degreeAllT[match(names(stat1),graphStats$genes)] =  stat1
  
  stat1 = graph.strength(gtumor,mode="out",loops=F)
  graphStats$degreeOutT[match(names(stat1),graphStats$genes)] =  stat1
  
  stat1 = graph.strength(gtumor,mode="in",loops=F)
  graphStats$degreeInT[match(names(stat1),graphStats$genes)] =  stat1
  
  print("btwn")
  btwCNT.t = betweenness(gtumor,weights = 1/E(gtumor)$weight)
  
  graphStats$betweenessT[match(names(btwCNT.t),graphStats$genes)] =  btwCNT.t
  
  print("eigen")
  eigIg.t = eigen_centrality(gtumor,weights = E(gtumor)$weight )$vector
  graphStats$eigIg.t[match(names(evcent(gtumor)$vector),graphStats$genes)] =  eigIg.t
  
  print("PR")
  #page rank
  pr.t = page_rank(gtumor,directed = F)
  graphStats$pageRankT[match(names(pr.t$vector),graphStats$genes)] =  pr.t$vector
  
  graphStats[is.na(graphStats)]=0
  
  
  write.table(graphStats,file=paste(folder,"/graph_stats_weight",label,".tsv",sep=""),quote = F,sep = "\t",row.names = F)
  
}


influ = function(aracneFile = "network.txt",genelist=NULL,cores=1){
  
    
  tumor = read.table(aracneFile,header = T,sep = "\t",stringsAsFactors = F)
  tumor$MI[tumor$ MI<=0]=1E-5
  if(!is.null(genelist)){
    tumor = tumor[which(tumor$Regulator %in% genelist & tumor$Target %in% genelist),]
  }
 
  gtumor = graph.data.frame(tumor)
  
  E(gtumor)$weight = E(gtumor)$MI
  
  
  
  source("collective_influence_algorithm.R")
  g = getInfluencers(gtumor,3,verbose = T,cores = cores)
  return(list(totalGenes = unique(c(tumor$Regulator,tumor$Target)),influencers=g$influencers))
}


metrics_weigthed_ARACNE2("../data/ARACNE_all_res/network.txt",folder = "../data/",label = "_tumor")
influencers = influ("../data/ARACNE_all_res/network.txt",cores=15)
save(influencers = file="../data/influencers.RData")

# to plot supfig 1 panel D
load("../data/complexes_v1.2_BRCA.RData")
genesInf = totalGenes %in% kk2$influencers
genesChrom = totalGenes %in% complexes$hgnc_symbol
fisher.test(table(genesChrom,genesInf))


# toplot supfig 1 panel c
library(ggplot2)

metrix = read.table("../data/graph_stats_weight_tumor.tsv",stringsAsFactors = F,header = T,sep = "\t")
metrix$complex = complexes$name[match(metrix$genes,complexes$hgnc_symbol)]
metrix.l = reshape2::melt(metrix)
metrix.l$complex = complexes$name[match(metrix.l$genes,complexes$hgnc_symbol)]


gois = intersect(complexes$hgnc_symbol,metrix$genes)
metrix$gois = metrix$genes %in% gois
library(ggpubr)
p1= ggdensity(metrix,x=c("degreeAllT"),combine=T,color="gois",fill="gois")+ scale_x_continuous(trans='log')+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG") )+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG"))+xlab("log degree")

p2 = ggdensity(metrix,x=c("betweenessT"),combine=T,color="gois",fill="gois")+ scale_x_continuous(trans='log')+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG"))+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG"))+xlab("log betweenness")

p3 = ggdensity(metrix,x=c("pageRankT"),combine=T,color="gois",fill="gois")+ scale_x_continuous(trans='log')+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG"))+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("nonCRG","CRG"))+xlab("log page rank")

ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggsave("../plots/SupFig1PanelC.pdf",width = 4,height = 9)


# toplot Supfig 1 panel b
evalCentr = function(gois,cent = "degreeAllT"){
  sum(metrix[metrix$genes %in% gois,cent],na.rm = T)
}
getP = function(s,s0,N){
  #https://blogs.sas.com/content/iml/2011/11/02/how-to-compute-p-values-for-a-bootstrap-distribution.html
  (1+sum(s >= s0))/(N+1)
}

chrom.degreeAll = evalCentr(gois,"degreeAllT")
randomGOIS = lapply(1:10000,function(x) sample(metrix$genes,length(gois)))
nullDist.degreeAll = unlist(lapply(randomGOIS,function(x) evalCentr(x,"degreeAllT")))
which(sort(c(nullDist.degreeAll,chrom.degreeAll),decreasing = T)==chrom.degreeAll)/10000
# 1
print("pval degree")
pval = getP(nullDist.degreeAll,chrom.degreeAll,10000)
#9.99E-5
print(pval)


id = order(c(nullDist.degreeAll,chrom.degreeAll,chrom.degreeAll),decreasing = T)
toPlot = data.frame(id = id,
                    value=c(nullDist.degreeAll,chrom.degreeAll,chrom.degreeAll)[id],
                    dist=c(rep("Null",length(nullDist.degreeAll)),"Chromatin","Chromatin")[id],stringsAsFactors = F)
toPlot$id = factor(toPlot$id,levels = toPlot$id)
p1=ggplot(toPlot,aes(x=value,fill=dist,color=dist))+
  #geom_bar(stat="identity")+
  #geom_histogram(aes(y=..density..),binwidth = 0.05)+
  geom_density(alpha=0.2,adjust=1,trim=T,aes(y=..scaled..))+xlab("Sum Degree")+
  theme_bw()+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))
  #scale_fill_discrete(name="",labels=c("CRG","Non CRG"))+scale_color_discrete(guide='none')



chrom.betweenessT = evalCentr(gois,"betweenessT")
randomGOIS = lapply(1:10000,function(x) sample(metrix$genes,length(gois)))
nullDist.betweenessT = unlist(lapply(randomGOIS,function(x) evalCentr(x,"betweenessT")))
which(sort(c(nullDist.betweenessT,chrom.betweenessT),decreasing = T)==chrom.betweenessT)
# 1
pval=getP(nullDist.betweenessT,chrom.betweenessT,10000)
print("pval betweeness")
print(pval)
#1.4E-3


id = order(c(nullDist.betweenessT,chrom.betweenessT,chrom.betweenessT),decreasing = T)
toPlot = data.frame(id = id,
                    value=c(nullDist.betweenessT,chrom.betweenessT,chrom.betweenessT)[id],
                    dist=c(rep("Null",length(nullDist.betweenessT)),"Chromatin","Chromatin")[id],stringsAsFactors = F)
toPlot$id = factor(toPlot$id,levels = toPlot$id)
p2= ggplot(toPlot,aes(x=value,fill=dist,color=dist))+scale_x_log10()+
  geom_density(alpha=0.2,adjust=1,trim=T,aes(y=..scaled..))+xlab("Sum Betweenness")+
  theme_bw()+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))
  

chrom.pageRankT = evalCentr(gois,"pageRankT")
randomGOIS = lapply(1:10000,function(x) sample(metrix$genes,length(gois)))
nullDist.pageRankT = unlist(lapply(randomGOIS,function(x) evalCentr(x,"pageRankT")))
which(sort(c(nullDist.pageRankT,chrom.pageRankT),decreasing = T)==chrom.pageRankT)
# 1
pval=getP(nullDist.pageRankT,chrom.pageRankT,10000)
print("pval page rank")
print(pval)
#1.E-4

id = order(c(nullDist.pageRankT,chrom.pageRankT,chrom.pageRankT),decreasing = T)
toPlot = data.frame(id = id,
                    value=c(nullDist.pageRankT,chrom.pageRankT,chrom.pageRankT)[id],
                    dist=c(rep("Null",length(nullDist.pageRankT)),"Chromatin","Chromatin")[id],stringsAsFactors = F)
toPlot$id = factor(toPlot$id,levels = toPlot$id)
p3=ggplot(toPlot,aes(x=value,fill=dist,color=dist))+scale_x_log10()+#scale_x_continuous(trans = scales::exp_trans())+
  geom_density(alpha=0.2,adjust=1,trim=T,aes(y=..scaled..))+xlab("Sum PageRank")+
  theme_bw()+
  scale_color_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))+
  scale_fill_manual(values=viridis::viridis_pal()(3)[-3],labels=c("CRG","Null"))

ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggsave("../plots/SupFig1PanelB.pdf",width = 4,height = 9)


