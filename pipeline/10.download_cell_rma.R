# download raw cel data
library(GEOquery)
ch.names= c("UPP","IRB","KAO","MAIRE","STK")
ch.geos = c("GSE3494","GSE45255","GSE20685","GSE65194","GSE1456")

tmpDir = "../data/clinical/cells"



for(i in 1:length(ch.geos)){
  print("********************************************************************")
  print(ch.geos[i])
  print("********************************************************************")
  expectedFilePrefixPattern="*"
  fileDataFrame = getGEOSuppFiles(ch.geos[i],makeDirectory=T, baseDir=tmpDir)    
  inFilePattern = ch.geos[i]
  if (substr(inFilePattern, 1, 3) == "GSE")
  {
    tarFilePath = file.path(tmpDir, paste(inFilePattern,"/",inFilePattern, "_RAW.tar", sep=""))
    
    if (!file.exists(tarFilePath))
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))
    
    downloadDir = paste(tmpDir,"/",inFilePattern,"/",sep="")
    dir.create(downloadDir, recursive=TRUE)
    untar(tarFilePath, exdir=downloadDir)
    inFilePattern = file.path(downloadDir, expectedFilePrefixPattern, sep="")
  }
  
  if (substr(inFilePattern, 1, 3) == "GSM")
  {
    downloadedFiles = list.files(path=tmpDir, full.names=TRUE, pattern=glob2rx(expectedFileSuffixPattern), ignore.case=TRUE)
    
    if (length(downloadedFiles) == 0)
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))
    
    inFilePattern = file.path(tmpDir, basename(downloadedFiles))
  }
  fileNamePattern = sub("\\-", "\\\\-", glob2rx(basename(inFilePattern)))
  fileNamePattern = sub("\\+", "\\\\+", fileNamePattern)
  celFilePaths = list.files(path=dirname(inFilePattern), pattern=fileNamePattern, full.names=TRUE, ignore.case=TRUE)
  celFilePaths = unique(celFilePaths)
  cellFilesPath2 = list.files(path=dirname(inFilePattern), pattern="CEL.gz", full.names=TRUE, ignore.case=TRUE)
  if(!ch.geos[i] %in% c("GSE3494")){
  celF = ReadAffy(filenames = fns)
  eset <- rma(celF)

  dir.create(paste("../data/clinical/rma/"),recursive = T)
  
  save(eset,file=paste("../data/clinical/rma/", inFilePattern,".RData",sep=""))
  }
  gc()
    
}




