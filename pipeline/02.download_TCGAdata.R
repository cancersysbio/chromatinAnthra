# Download TCGA breast cancer data gdc.cancer.gov

library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      file.type="counts.gz")
    
GDCdownload(query, method = "api", files.per.chunk = 2,directory = "../data")
data.exp <- GDCprepare(query,directory = "../data")
save(data.exp,file=paste("../data/02.GeneExp_gene_counts.RData",sep=""))

query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      file.type="FPKM.txt.gz")
    
GDCdownload(query, method = "api", files.per.chunk = 2,directory = "../data")
data.exp <- GDCprepare(query,directory = "../data",remove.files.prepared=F)
save(data.exp,file="../data/02.GeneExp_gene_FPKM.RData")

query <- GDCquery(project = projects.TCGA[i],
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    file.type="FPKM-UQ.txt.gz")

GDCdownload(query, method = "api", files.per.chunk = 2,directory = "../data")
data.exp <- GDCprepare(query,directory = "../data",remove.files.prepared=F)
save(data.exp,file="../data/02.GeneExp_gene_FPKM-UQ.RData")