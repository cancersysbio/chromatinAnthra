# generate Chromatin Regulator Genes list

# complexes file (v1.2) was generated on Feb 2018. Changes in GO modify the final list
# of CRGs. For reproductibility purposes, please use file 1.2

library(biomaRt)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


complexes = NULL

# methyl transferases
HMT <- getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0018024', mart = ensembl)
HMT = HMT[HMT$go_id=='GO:0018024',]
HMT = HMT[HMT$chromosome_name %in% c(1:22,"X","Y"),]

HMT$name="METHTRANSF"
HMT$color="#D9D9D9"

complexes = rbind(complexes,HMT)

# histone demethylases
HDMT = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0032452', mart = ensembl)
HDMT = HDMT[HDMT$go_id=='GO:0032452',]
HDMT = HDMT[HDMT$chromosome_name %in% c(1:22,"X","Y"),]

HDMT$name="DEMETH"
HDMT$color="#FFED6F"

complexes = rbind(complexes,HDMT)

# histone deacetylase
HDA =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
              filters = 'go', values = 'GO:0004407', mart = ensembl)
HDA = HDA[HDA$go_id=='GO:0004407',]
HDA = HDA[HDA$chromosome_name %in% c(1:22,"X","Y"),]

HDA$name="DEACETH"
HDA$color="darkorchid1"

complexes = rbind(complexes,HDA)

# histone acetyltransferase
HA =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0004402', mart = ensembl)
HA = HA[HA$go_id=='GO:0004402',]
HA = HA[HA$chromosome_name %in% c(1:22,"X","Y"),]

HA$name="ACETH"
HA$color="cornsilk"

complexes = rbind(complexes,HA)

# histone phosporilation
HP =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0016572', mart = ensembl)
HP = HP[HP$go_id=='GO:0016572',]
HP = HP[HP$chromosome_name %in% c(1:22,"X","Y"),]

HP$name="PHOSPH"
HP$color="brown1"

complexes = rbind(complexes,HP)

# PRC1 comples
PRC1 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0035102', mart = ensembl)
PRC1 = PRC1[PRC1$go_id=='GO:0035102',]
PRC1 = PRC1[PRC1$chromosome_name %in% c(1:22,"X","Y"),]
PRC1$name="PRC1"
PRC1$color="#FFFFB3"

complexes = rbind(complexes,PRC1)

# PRC2 complex
PRC2 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0035098', mart = ensembl)
PRC2 = PRC2[PRC2$go_id=='GO:0035098',]
PRC2 = PRC2[PRC2$chromosome_name %in% c(1:22,"X","Y"),]

PRC2$name="PRC2"
PRC2$color="#BEBADA"

complexes = rbind(complexes,PRC2)

# REMODELERS FAMILY according to http://www.nature.com/nrm/journal/v7/n6/full/nrm1945.html
#SWI/SNF, ISWI, NURD/Mi-2/CHD, INO80 and SWR1.

# BAF
# Tumour suppressor, Differentiation, Development, Elongation, Signalling, Splicing

BAF =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
              filters = 'go', values = 'GO:0016514', mart = ensembl)
BAF = BAF[BAF$go_id=='GO:0016514',]
otherBAF = c("PBRM1","SS18","BCL11A","BCL11B","BCL7A","BCL7B","BCL7C","DPF1","DPF2","DPF3","PHF10","BRD9")
BAF2 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = otherBAF, mart = ensembl)
BAF2$go_id = BAF$go_id[1]
BAF2$name_1006 = BAF$name_1006[1]
BAF=rbind(BAF,BAF2)
BAF = BAF[BAF$chromosome_name %in% c(1:22,"X","Y"),]

BAF$name = "SWI_SNF"
BAF$color = "#8DD3C7"

complexes = rbind(complexes,BAF)

## ISWI  (NURF / ACF / CHRAC  /WICH / NORC /RSF / CERF)
# NURF: SMARCA1 + BPTF
# WICH: SMARCA5 + WSTF
# NORC: SMARCA5 + TIP5 (BAZ2A) + BAZ2B
# RSF: SMARCA5 + RSF1
# ACF: ACF1 + SMARCA5
# CHRAC: SMARCA5 + ACF1 + CHRAC15 + CHRAC17 
# CERF: SMARCA1 + CECR2

iswi_genes = c("SMARCA5","SMARCA1","BPTF","BAZ1B","BAZ2A","BAZ2B","RSF1","BAZ1A","CHRAC1","POLE3","CECR2","C17orf49","RBBP4","HMGXB4","USF1")
ISWI = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = iswi_genes, mart = ensembl)

ISWI$go_id = ""
ISWI$name_1006 = "ISWI"
ISWI = ISWI[ISWI$chromosome_name %in% c(1:22,"X","Y"),]

ISWI$name="ISWI"
ISWI$color = "#FDB462"

complexes = rbind(complexes,ISWI)

# CHDs NURD/Mi-2
# Transcriptional repression and silencing, Development

chdsGenes = c("CHD1","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9","HDAC1","HDAC2","MBD2","MBD3","MTA1", "MTA2","MTA3","GATAD2A","GATAD2B", "RBBP4","RBBP7")
CHD = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
            filters = 'hgnc_symbol', values = chdsGenes, mart = ensembl)

CHD$go_id = ""
CHD$name_1006 = "CHD_NURD_mi2"
CHD = CHD[CHD$chromosome_name %in% c(1:22,"X","Y"),]

CHD$name="CHD_NURD"
CHD$color = "#FB8072"

complexes = rbind(complexes,CHD)

# INO80
#
# http://www.nature.com/nrm/journal/v10/n6/fig_tab/nrm2693_T1.html
INO80 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
                filters = 'go', values = 'GO:0031011', mart = ensembl)
INO80 = INO80[INO80$go_id=='GO:0031011',]

ino80_genes = c("ACTB")
INO80_2 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = ino80_genes, mart = ensembl)
INO80_2$go_id = INO80$go_id[1]
INO80_2$name_1006 = INO80$name_1006[1]
INO80=rbind(INO80,INO80_2)
INO80 = INO80[INO80$chromosome_name %in% c(1:22,"X","Y"),]

INO80$name="INO80"
INO80$color="cadetblue2"

complexes = rbind(complexes,INO80)

# SWR1
# DNA repair
# http://www.nature.com/nrm/journal/v10/n6/fig_tab/nrm2693_T1.html
swr_genes = c("SRCAP","RUVBL1","RUVBL2","ACTB","ACTL6A","ACTR6","YEATS4","DMAP1","ERCC5","VPS72")
SWR1 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = swr_genes, mart = ensembl)
SWR1 = SWR1[SWR1$chromosome_name %in% c(1:22,"X","Y"),]

SWR1$go_id = ""
SWR1$name_1006 = "SWR1 complex"

SWR1$name="SWR1"
SWR1$color="goldenrod4"

complexes = rbind(complexes,SWR1)

# BAP1/PR-DUB
bap1_bpdub_complex_genes = c("ASXL1","ASXL2","BAP1","BARD1","BRCA1","FOXK1","FOXK2","HCFC1","KDM1B","OGT","YY1")
BAP1_PRDUB = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                   filters = 'hgnc_symbol', values = bap1_bpdub_complex_genes, mart = ensembl)
BAP1_PRDUB$go_id=""
BAP1_PRDUB$name_1006="BAP1_PR-DUB_Complex"
BAP1_PRDUB = BAP1_PRDUB[BAP1_PRDUB$chromosome_name %in% c(1:22,"X","Y"),]

BAP1_PRDUB$name="BAP1_PRDUB"
BAP1_PRDUB$color="#B3DE69"

complexes = rbind(complexes,BAP1_PRDUB)


# CAF-1
CAF1 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0033186', mart = ensembl)

CAF1 = CAF1[CAF1$go_id=='GO:0033186',]

CAF1$name="CAF1"
CAF1$color="#80B1D3"

complexes = rbind(complexes,CAF1)

# cohesin 
# http://www.nature.com/nrc/journal/v14/n6/full/nrc3743.html
#http://www.nature.com/nrg/journal/v15/n4/full/nrg3663.html
# check role CTCF in  methylation
cohe_genes = c("STAG1","STAG2","STAG3","SMC1A","SMC1B","SMC3","RAD21","REC8","RAD21L1","WAPL","PDS5A","PDS5B","CDCA5","MTA3")
COHE =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
               filters = 'hgnc_symbol', values = cohe_genes, mart = ensembl)
COHE$go_id=""
COHE$name_1006="Cohesin_Complex"

COHE$name="COHE"
COHE$color="#BC80BD"

complexes = rbind(complexes,COHE)

# condensins
#http://genesdev.cshlp.org/content/26/15/1659.full
conde_genes = c("SMC4","NCAPH","NCAPG","NCAPD2","SMC2","NCAPD3","NCAPG2","NCAPH2")
CONDE =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = conde_genes, mart = ensembl)
CONDE$go_id=""
CONDE$name_1006="Condensin_Complex"

CONDE$name="CONDE"
CONDE$color="#CCEBC5"

complexes = rbind(complexes,CONDE)

#TOPo activity
TOPO = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0003916', mart = ensembl)
TOPO = TOPO[TOPO$go_id=='GO:0003916',]

TOPO$name="TOPO"
TOPO$color="#FCCDE5"

complexes = rbind(complexes,TOPO)


#  DNA methyltransferases

DNA_meth = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
                   filters = 'go', values = 'GO:0006306', mart = ensembl)
DNA_meth = DNA_meth[DNA_meth$go_id=='GO:0006306',]
DNA_meth$name="DNA_METH"
DNA_meth$color="aquamarine"

complexes = rbind(complexes,DNA_meth)


#  DNA demethylases
DNA_demeth = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
      filters = 'go', values = 'GO:0080111', mart = ensembl)
DNA_demeth = DNA_demeth[DNA_demeth$go_id=='GO:0080111',]

DNA_demeth$name="DNA_DMETH"
DNA_demeth$color="slateblue1"

complexes = rbind(complexes,DNA_demeth)

hist_genes = read.table("../data/histone_genenames.txt",header = T,sep = "\t",stringsAsFactors = F)

HIST =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
               filters = 'hgnc_symbol', values = hist_genes$Approved.Symbol, mart = ensembl)
HIST$go_id=""
HIST$name_1006="HISTONES"

HIST$name="HIST"
HIST$color="tan3"
complexes = rbind(complexes,HIST)


pio_genes = c("FOXA1")

pio =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = pio_genes, mart = ensembl)
pio$go_id=""
pio$name_1006="pio_GENES"

pio$name="pio"
pio$color="violetred1"

complexes = rbind(complexes,pio)


# remove duplicate entrez
complexes2 = unique(complexes[order(complexes$entrezgene_id),c(1,2)])
idDup = duplicated(complexes2$hgnc_symbol)
falseEntrez = complexes2$entrezgene_id[idDup]

complexes = complexes[!complexes$entrezgene_id %in% falseEntrez,]
complexes = na.omit(complexes) # 471 unique genes



save(complexes, file="../data/complexes_v1.3_BRCA.RData")
