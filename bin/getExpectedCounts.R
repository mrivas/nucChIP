rm(list=ls())
set.seed(1)
library(scales)
library(fields)
library(graphics)
library(doParallel)
cores=makeForkCluster(3)
registerDoParallel(cores)
print("hello")
getData = function(file){
	classes <- rep("NULL",11)
	classes[9]="numeric" # -log(p-value,10)
	classes[11]="numeric" # counts
	data=read.table(file,skip = 1, colClasses = classes)
	counts=subset(data[[2]],data[[1]]>=0)
	# mark outliers
	threshold=quantile(counts,0.999)
	counts[counts>threshold] = NA
	#counts=squish(counts,quantile(counts,c(0,.99))  ) # cap outliers
	return(counts)
}

print("definefunction")
expectedValues <- function(file,mnase){

	print( file )
	# Load data
	name=unlist(strsplit(tail(unlist(strsplit(file,"/")),1) , "\\."))[1]
	data=getData(file)
	
	fudgeit <- function(){
	  xm <- get('xm', envir = parent.frame(1))
	  ym <- get('ym', envir = parent.frame(1))
	  #z  <- get('dens', envir = parent.frame(1))
	  colramp <- get('colramp', parent.frame(1))
	  numdup=aggregate(numdup~.,data.frame(mnase,data,numdup=1),length)
	  image.plot(zlim=c(0,max(numdup$numdup)), col = colramp(256), legend.only = T, add =F)
	}
	
	# Aggregate data
	agdata=aggregate(data,by=list(mnase),FUN=mean,na.rm=T)
	names(agdata)=c("MNase_counts",paste(name,"_expected_counts",sep=""))
	# Linear transformation
	reference = agdata[1]*sum(data,na.rm=T)/sum(mnase,na.rm=T)
	# Save image
	ag_lower=aggregate(data,by=list(mnase),FUN=quantile,probs=0,na.rm=T)
	ag_upper=aggregate(data,by=list(mnase),FUN=quantile,probs=0.95,na.rm=T)
	lower=ag_lower[[2]]
	upper=ag_upper[[2]]
	svg(paste(name,".expectedCounts.svg",sep=""),width = 5, height = 7)
	par(mar = c(5,4,4,5) + .1)
	smoothScatter(mnase,data,xlab="MNase counts per nucleosome",ylab=paste("Histone counts per nucleosome"),main=name,nrpoints=0,postPlotHook=fudgeit)
	lines(agdata,col="blue")
	lines(upper,col="red")
	lines(lower,col="red")
	lines(agdata[[1]],reference[[1]],col="black",lty=2)
	legend("topleft", c("Expected histone counts", "95% CI","Linear transformation"), col=c("blue","red","black"),lty=c(1,1,2),bty="n" )
	dev.off()
	# Save aggregated data as tables
	newData=cbind(agdata,upper)
	fName=paste(name,".expectedCounts.txt",sep="")
	write.table(newData,file=fName,quote=FALSE,sep="\t",row.names=F,col.names=T)
}

print("histogra")
# Mnase counts
mnase=getData("/data2/rivasas2/singleNucleosome/secondBatch/counts/allData/8_mnase.allDataNuc.counts.bed")
svg("8_mnase.svg")
hist(mnase,xlab="MNase counts per nucleosome",main="8_mnase")
dev.off()
# MNase ChIP-seq counts
print("obtainfiles")
folder="/data2/rivasas2/singleNucleosome/secondBatch/combineLib/allData/"
#pattern="*.8_mnase.counts.bed"
pattern="*.allDataNuc.counts.bed"
files=c()
for( file in list.files(path=folder,pattern=pattern)  ) {
	if( substr(file,1,7)!="8_mnase"){
		files=c(files,paste(folder,file,sep=""))
	}
}
print("Running in parallel")
foreach( i=1:length(files) ) %dopar% expectedValues(files[i],mnase)
