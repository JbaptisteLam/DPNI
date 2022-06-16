##################################################################################################################################################
##################################################################################################################################################
##### SeqFF.R
##### Sung Kim, Phd
##### Sequenom, Inc.
##### 
##################################################################################################################################################
##### Instructions
##### Command Line 
##### R CMD BATCH 
##### --i input directory 
##### --f input file name
##### --d output directory
##### --o output file name
##### --t data type; sam file (without header) or tabulated read counts ordered by genomic coordinates found in SupplementalTable1.csv
##### SeqFF.R


ff.pred <- function(gc.norm.bc.61927, B, mu, parameter.1, parameter.2){
  gc.norm.bc.61927[is.na(gc.norm.bc.61927)] <- 0
  gc.norm.bc <- gc.norm.bc.61927[grepl('chr[0-9]', names(gc.norm.bc.61927))]
  gc.norm.bc.c <- gc.norm.bc - mu
  y.hat <- matrix(c(1, gc.norm.bc.c), nrow = 1) %*% B
  y.hat.rep <- sum(y.hat, na.rm = T) / sum(gc.norm.bc)
  ff.hat <- (y.hat.rep + parameter.1) * parameter.2
  return(ff.hat)
  
}

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##### COMMAND LINE ARGUMENTS
arg=commandArgs()
input.dir = unlist(strsplit(arg[ pmatch("--i",arg)], "="))[2]
file.name = unlist(strsplit(arg[ pmatch("--j",arg)], "="))[2]
output.dir = unlist(strsplit(arg[ pmatch("--d",arg)], "="))[2]
output.filename = unlist(strsplit(arg[ pmatch("--o",arg)], "="))[2]
datatype = unlist(strsplit(arg[ pmatch("--t",arg)], "="))[2]


load("SupplementalFile1.RData")
bininfo = read.csv("SupplementalTable2.csv")
colnames(bininfo)[1]="binName"
bininfo$binorder=c(1:61927)

setwd(input.dir)
print(input.dir)
print(file.name)
print(output.dir)
print(output.filename)
print(datatype)
##################################################################################################################################################
##### READ IN DATA
if( datatype =="sam" )
{

dat <- read.table(pipe(paste("cut -f3,4",file.name,sep=" ")), header=FALSE,colClasses= c("character","integer"))
colnames(dat)=c("refChr","begin")


dat=dat[dat$refChr!="*" & dat$refChr!="chrM" ,]
print(head(dat))

#modulo et trunc
binindex = ifelse((dat$begin%%50000)==0,floor(dat$begin/50000)-1,floor(dat$begin/50000))
print(head(binindex))

fastbincount = table(paste(dat$refChr,binindex, sep="_"))
print(fastbincount)
newtemp=as.data.frame(matrix(NA, ncol=2, nrow=length(fastbincount)))
newtemp[,1]=paste(names(fastbincount))
newtemp[,2]=as.numeric(paste(fastbincount))
colnames(newtemp)=c("binName","counts")

}
if( datatype =="counts" )
{
newtemp <- read.table(file.name, header=FALSE,colClasses= c("character","integer"))
colnames(newtemp)=c("binName","counts")
}

bininfo = merge(bininfo, newtemp, by="binName",all.x=T)
bininfo=bininfo[order(bininfo$binorder),]
print(head(newtemp))
print(head(bininfo))
##################################################################################################################################################
##### DATA PROCESSING
autosomebinsonly = bininfo$BinFilterFlag==1 & bininfo$CHR!="chrX" & bininfo$CHR!="chrY"
alluseablebins = bininfo$BinFilterFlag==1 
autoscaledtemp  <- bininfo$counts[autosomebinsonly]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
allscaledtemp  <- bininfo$counts[alluseablebins]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
# additive loess correction
mediancountpergc <- tapply(
autoscaledtemp,bininfo$GC[autosomebinsonly], function(x) median(x, na.rm=T))
## prediction 
loess.fitted  <- predict( loess(mediancountpergc ~ as.numeric(names(mediancountpergc))), bininfo$GC[alluseablebins]) 
normalizedbincount <- allscaledtemp + ( median(autoscaledtemp, na.rm=T) - loess.fitted )

bincounts=rep(0,61927)
names(bincounts) = bininfo$binName
bincounts[alluseablebins] <- (normalizedbincount/sum(normalizedbincount, na.rm=T)) * length(normalizedbincount)

wrsc=ff.pred(bincounts,B,mu,parameter.1,parameter.2)
enet = bincounts %*% elnetbeta+elnetintercept
ff=c((wrsc+enet)/2, enet, wrsc)
names(ff)=c("SeqFF","Enet","WRSC")

setwd(output.dir)

write.csv(ff, file=output.filename)



