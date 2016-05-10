#Set working directory to file directory with data
setwd("E:/Noble Foundation/Output3")

#Set the desired volume cutoff (cubic microns)
minVolumeCutoff<-2

#This function counts the total number of bundles in each file and returns the count
countObjects<-function(fileName){
	#Load the file into a temporary matrix
	fileMatrix<-read.table(fileName, sep="\t", header=TRUE)
	
	#Return the number of samples in the file
	count<-nrow(fileMatrix)
	return(count)
}

#This function fills a matrix with the corresponding measurements in an excel file
excelToMatrix<-function(fileName, maxCount){
	#Load the file into a temporary matrix
	fileMatrix<-read.table(fileName, sep="\t", header=TRUE)
	
	#Retrieve the volume measurements
	volumeVector<-fileMatrix$Volume
	
	#Fill the remaining vector with NA so that sapply creates a matrix (all returned vectors are of equal length)
	for (a in length(volumeVector):maxCount){
		volumeVector[a+1]<-NA
	}
	
	#Return the volume measurements column
	return(volumeVector)
}

#Get all file names that end with a *.xls file extension
fileNameList<-list.files(pattern = "\\.xls$",  ignore.case=TRUE)

#Count the number of results in each file
fileObjectCounts<-sapply(fileNameList, countObjects, simplify = TRUE)

#Find the largest number of object counts (this value will be used to setup a matrix to store all volume measurements)
maxCount<-max(fileObjectCounts)

#Add the volume measurements to  each corresponding column in the a new matrix
allVolumeMatrix<-sapply(fileNameList, excelToMatrix, maxCount, simplify = TRUE, USE.NAMES = FALSE)

#Extract the sample ID for each file by removing the prefix and extension and set as the column names in the new matrix
sampleIDvector<-sapply(fileNameList, gsub, pattern="Statistics for ", replacement="")
sampleIDvector<-sapply(sampleIDvector, gsub, pattern=".xls", replacement="")
colnames(allVolumeMatrix)<-sampleIDvector

#Remove the first column since it is the summary file
allVolumeMatrix<-allVolumeMatrix[,-1]

#Create a vector of all sample IDs so that replicates can be pooled together
#initize variables and vectors
index<-1 
genotypeNameVector<-vector(,0)
genotypeVolumeList<-list()

#Remove the fist item from the name vector since it is no longer used
sampleIDvector<-sampleIDvector[-1]

for(a in 1:length(sampleIDvector)){
	#Remove the numeric sample identifier suffix from the name
	genotype<-gsub(" - [0-9]+", "", sampleIDvector[a])
	#On the first iteration, add the genotype and count up the index
	if(a == 1){
		genotypeNameVector[index]<-genotype
		genotypeVolumeList<-allVolumeMatrix[,a]
	}
	
	#On all following iterations, check to see if the same genotype is found, and if so, add it to the vector 
	else if(genotype == genotypeNameVector[index]){
		genotypeVolumeList<-c(genotypeVolumeList, allVolumeMatrix[,a])
	}
	
	#If a new genotype is found, create a new list and add it to the list of lists
	else{
		#If index is 1 then create a new list of lists
		if(index == 1){
			listOfLists<-list(genotypeVolumeList)
		}
		#Otherwise, add the current list to the list of lists
		else{
			listOfLists[[length(listOfLists)+1]]<-genotypeVolumeList
		}
		
		#Increase the index counter
		index<-index+1
		genotypeNameVector[index]<-genotype
		genotypeVolumeList<-allVolumeMatrix[,a]
				
	}
}
#Add the last list to the list of lists
listOfLists[[length(listOfLists)+1]]<-genotypeVolumeList

#convert the list of lists to a matrix
genotypeVolumeMatrix<-do.call(cbind, lapply(listOfLists, unlist))
colnames(genotypeVolumeMatrix)<-genotypeNameVector

#Remove all volume measurements that are below the set cutoff
genotypeVolumeMatrix<-replace(genotypeVolumeMatrix, genotypeVolumeMatrix < minVolumeCutoff, NA)

#Create a notched box plot for each genotype
par(mar = c(15,6,3,3))
boxplot(genotypeVolumeMatrix, log="y", las = 2, names = genotypeNameVector, notch = TRUE, ylim = c(min(genotypeVolumeMatrix, na.rm = TRUE), max(genotypeVolumeMatrix, na.rm = TRUE)), ps = 1, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1, bty="n")
mtext(expression(paste( plain("GFP Volume Distribution (μm") ^ plain("3"), plain(")") )), side=2, line = 3, cex = 1)
