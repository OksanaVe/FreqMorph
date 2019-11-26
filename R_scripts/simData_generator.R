library(diversitree)
library(phangorn)
library(greenbrown)

# Get full paths for all True trees in directory
trees <- list.files("/FreqMorph/SimData/True_trees", pattern="*.tre", full.names=TRUE)

setwd("/FreqMorph/SimData")


# Use Liam Revel's wrapper to write matrices in standard data format; modified to specify symbols for compatibility with software other than Mesquite
writeNexusData<-function(x, file, format = "dna", interleaved=FALSE){
  if(format=="dna"||format=="DNA") write.nexus.data(x,file,format,interleaved=FALSE)
  else if(format=="protein"||format=="PROTEIN") write.nexus.data(x,file,format,interleaved=FALSE)
  else if(format=="standard"||format=="STANDARD"){
   X<-vector(mode="list",length=nrow(x))
   for(i in 1:nrow(x)) X[[i]]<-x[i,]
   names(X)<-rownames(x)
   write.nexus.data(X,file,format="dna",interleaved=FALSE)
   ff<-readLines(file)
   ff<-gsub("DNA","STANDARD SYMBOLS=\"1 2\" ",ff)
   write(ff,file)
  }
}

# Read each tree and generate data matrix
for (n in 1:100) {                     # this sets number of data set replicates per tree
for (i in 1:length(trees))            # create data matrix per each true tree
{ tre <- read.tree(trees[i])          # read each true tree (iteration)
chars <- 20*Ntip(tre)                # set number of characters in data matrix to 2.5 times number of taxa
spec <- paste0(basename(trees[i]), "_", n, ".nex")  #define data matrix file name (Tree_replicate.nex)
data=matrix(0,Ntip(tre),chars)       # create emptry matrix with specified dimensions
rownames(data) <- tre$tip.label      # add taxa names to the matrix (matching to the tree tip labels)
for (j in 1:ncol(data))  {           # fill in data matrix
p=runif(1)
q=1-p
lambda=runif(1, min=5, max=500)
char <- rTraitDisc(tre, model="ARD", k=2, rate=rexp(2, lambda), freq=c(p, q)) # generate each character

	if(AllEqual(char) == TRUE) {
	repeat {
p=runif(1)
q=1-p
lambda=runif(1, min=5, max=500)
char <- rTraitDisc(tre, model="ARD", k=2, rate=rexp(2, lambda), freq=c(p, q))
data[, j]=char
if (AllEqual(char) == FALSE) { break}
}	
					} else {data[, j]=char}
}



dataPhy <- phyDat (data, type="USER", levels=c("0", "1", "2", "3"))
cindex <- CI(tre, dataPhy, sitewise=FALSE)
cindex
writeNexusData(data, spec, format="standard", interleaved=FALSE) # write data matrix to file before next iteration
write.table(c(spec, cindex), "CI_data.txt", append=TRUE)
}
}

