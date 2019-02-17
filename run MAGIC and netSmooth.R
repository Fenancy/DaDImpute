require(SingleCellExperiment)
require(netSmooth)
library(Rmagic)
library(viridis)
library(phateR)

############initialisation############
# number of iterations to run
n<-1
#read scRNA seq (gene*cell matrix)
load("../subraw.rda")
raw<-subraw
#trim off the version number in gene ensembl code if needed
rownames(raw)<-strtrim(rownames(raw),15)
#load appropriate ppi
load("./human.ppi.rda")
ppi<-human.ppi
#get total number of non zeros
non0<-sum(raw!=0)
#get the indices .of non-zero entries
ind_non0<-which(raw!=0)
# p% non-zeros to be hidden 
for (p in c(99)){
	#calculate x: the number of entries to be masked
	x<-as.integer(non0*p/100)
	#matrix to store Spearman rank coeffecients
	rho<-matrix(nrow=2,ncol=n,dimnames=list(c("netSmooth","MAGIC"),1:n))
	#lists to store counts to be masked and imputed masked values (raw and logged)
	mask_mg<-mask_net<-mask<-list()
	for (i in (1:n)){
		#duplicate data
		dup<-raw
		#subsample the non-zero indices to be masked as zeros
		set.seed(i);ind_mask<-sample(ind_non0,size=x)
		#before masking, save the original counts
		mask[[i]]<-mask_raw<-raw[ind_mask]
		#masking
		dup[ind_mask]<-0
	 	#######################netSmoothing################
	 	#convert to singleCellExperiment object
		seq<-SingleCellExperiment(assays = list(counts = dup))
		seq.sm.se <-netSmooth(seq, ppi, alpha=0.5)
		seq.sm.sce <-SingleCellExperiment(assays=list(counts=assay(seq.sm.se)), colData=colData(seq.sm.se))
		seq_net<-assay(seq.sm.sce)
		mask_net[[i]]<-seq<-seq_net[ind_mask]
		#Spearman correlation
	 	###########################Magic####################
	 	#tranpose the matrix so MAGIC accepts it
		tr<-t(dup)
		#Normalising Data
		tr<-library.size.normalize(tr)
		tr<-sqrt(tr)
		#Running Magic#
		tr <- magic(tr, genes="all_genes")
		tr <- data.matrix(tr$result)
	 	#tranpose the matrix back
	 	seq_mg<-t(tr)
	 	mask_mg[[i]]<-seq<-seq_mg[ind_mask]
	 	#Spearman correlation
		#######save imputed data 
		#fn_full=paste("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/CV analysis/full_",p,"percent.rda",sep="")
		#save(dup,seq_net,seq_mg,file=fn_full)
		}
	rho_net<-mean(rho[1,])
	rho_mg<-mean(rho[2,])
	fn_rho=paste("./rho_",p,"percent.rda",sep="")
	save(rho,rho_net,rho_mg,file=fn_rho)
	fn=paste("./missing_",p,"percent.rda",sep="")
	save(mask_mg,mask_net,mask,file=fn)
	fn_log=paste("./missing_",p,"percent_logged.rda",sep="")
	save(mask_mg_log,mask_net_log,mask_raw_log,file=fn_log)
	}

	