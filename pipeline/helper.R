# helper functions

# merge samples without batch effect removal
mergeNONE = function(esets)
{
	eset1 = esets[[1]];
	annot1 = annotation(eset1)
	
	for(i in 2:length(esets))
	{
		eset2 = esets[[i]];
		d1 = exprs(eset1);
		d2 = exprs(eset2);
		
		#-----------------------------------------------------
		# Rebuild fData
		#-----------------------------------------------------
		cg = sort(intersect(rownames(d1), rownames(d2))); 
		# If too few overlapping genes
		if(length(cg) < (min(dim(d1)[1],dim(d2)[1])/100))
		{
			msg(" ! WARNING ! Number of common genes < 1%")
		}
		fData = fData(eset1)[cg,]
		
		#-----------------------------------------------------
		# Rebuild pData
		#-----------------------------------------------------
		
		# Create new pData: first common annotations, then unique ones
		p1 = pData(eset1)
		p2 = pData(eset2)
		cp = sort(intersect(colnames(p1),colnames(p2)));
		tp = sort(unique(union(colnames(p1),colnames(p2))))
		
		sp1 = setdiff(colnames(p1),cp)
		sp2 = setdiff(colnames(p2),cp)
		
		pheno = matrix(NA,ncol=length(tp),nrow=nrow(p1)+nrow(p2))
		rownames(pheno) = c(rownames(p1),rownames(p2))
		colnames(pheno) = tp;
		
		if(length(cp)!=0)
		{
			pheno[1:nrow(p1),cp] = as.matrix(p1[,cp])
			pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)),cp] = as.matrix(p2[,cp])
		}
		
		if(length(sp1)!=0)
		{
			pheno[1:nrow(p1),sp1] = as.matrix(p1[,sp1])
		}
		if(length(sp2)!=0)
		{
			pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)), sp2] = as.matrix(p2[,sp2])
		}
		
		pData = as.data.frame(pheno,stringAsFactors=F)
		
		#-----------------------------------------------------
		# Rebuild rest eset
		#-----------------------------------------------------
		
		# drop = FALSE prevents to loose the column name in case there is only 1 sample in one exprs
		d1 = d1[cg, , drop = FALSE];
		d2 = d2[cg, , drop = FALSE];
		
		
		#eset1 = new("ExpressionSet");    
		eset1 = ExpressionSet(assayData =cbind(d1, d2) )
		#exprs(eset1) = cbind(d1, d2);
		pData(eset1) = pData;
		fData(eset1) = fData;
		
		annot1 = c(annot1, annotation(eset2));
		
	}
	
	annotation(eset1) = unique(annot1);
	return(eset1);
}

removeBatchCombat = function(dataset,batch){
  
  idna = which(is.na(batch))
  
  if(length(idna)>0){
    dataset2 = dataset[,-idna]
    batch2 = batch[-idna]
  }else{
    dataset2 = dataset
    batch2= batch
  }
  
  exp.data = exprs(dataset2)
  
  #mod = model.matrix(~as.factor(outcome[-idna]))
  mod0 = model.matrix(~1,data=data.frame(batch2,stringsAsFactors = F))
  
  modcombat = model.matrix(~1, data=data.frame(batch2,stringsAsFactors = F))
  combat_edata = sva::ComBat(dat=exp.data, batch=batch2, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  #pValuesComBat =  sva::f.pvalue(combat_edata,mod,mod0)
  #qValuesComBat = p.adjust(pValuesComBat,method="BH")
  
  
  #plot(-log10(qValuesComBat))
  #sum(qValuesComBat<0.05)
  
  dataset3 = dataset2
  exprs(dataset3)=combat_edata
  
  return(dataset3)
}