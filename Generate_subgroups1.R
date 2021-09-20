setwd("/nfs/proj/REPO-TRIAL/montserrat")

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Extract the optimal number of clusters (patient subgroups) per disease based on differential gene expression information ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration:
library("cluster");library("NbClust");library(e1071)
set.seed(1)
extract_subgroups<-function(matinter,enfermedad){
  k.max <- length(rownames(matinter))-1  ## Number of patients -1
  if(k.max>1){
    data <- matinter
    sil <- rep(0, k.max)
    # Compute the average silhouette width for
    # k = 2 to k = 15
    for(i in 2:k.max){
      tryCatch({km.res <-try( kmeans(data, centers = i, nstart = 100))
      ss <- silhouette(km.res$cluster, dist(data))
      sil[i] <- mean(ss[, 3])
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    # Plot the  average silhouette width
    plot(1:k.max, sil, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters k",main=enfermedad)
    abline(v = which.max(sil), lty = 2)
    pdf(file=paste("Data/Genes_based_clustering/Silhouette/",enfermedad,".pdf",sep=""))
    plot(1:k.max, sil, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters k",main=enfermedad)
    abline(v = which.max(sil), lty = 2)
    dev.off()
    km <- kmeans(matinter, centers = which.max(sil), nstart = 100)
    return(km)
  }
}

## We provide intervals (lines 49-53) for running the script from lines 36 to 67 in 5 cores at the same time
# args=commandArgs(trailingOnly = TRUE)
# argumento<-as.numeric(args)
# print(class(argumento))

starting<-Sys.time()
load("Data/Matrix_with_patients_and_genes_1_0_-1.Rdata") #newmatrix
# diseases<-gsub("_GSE.+","",rownames(newmatrix))
# disease<-unique(diseases)
info<-read.csv2("Patients_information_table_v2.txt",stringsAsFactors = F,sep=",",header=T)
diseases<-info[,2]
numbers<-table(info[,2])
numbers<-numbers[order(numbers,decreasing = F)]
disease<-names(numbers)


if("Genes_based_clustering"%in%list.files("Data/")==FALSE){dir.create("Data/Genes_based_clustering")}
if("Silhouette"%in%list.files("Data/Genes_based_clustering")==FALSE){dir.create("Data/Genes_based_clustering/Silhouette")}
if("Clusters"%in%list.files("Data/Genes_based_clustering")==FALSE){dir.create("Data/Genes_based_clustering/Clusters")}

# Remove from the loop the diseases already analysed
dones <- gsub(".pdf","", list.files("Data/Genes_based_clustering/Silhouette"))
if(length(dones) > 0) {
  disease <- setdiff(disease,dones)
}


print("We start the loop!")
for(a in 1:length(disease)){
  star<-Sys.time()
  dis<-disease[a]
  matinter<-newmatrix[grep(paste("^",dis,"$",sep=""),diseases),]
  ends<-Sys.time()
  clusters<-extract_subgroups(matinter,dis)
  save(clusters,file=paste("Data/Genes_based_clustering/Clusters/",dis,".Rdata",sep=""))
  print(paste(dis," finished, ",a," out of ",length(disease)," in ",ends-star,sep=" "))
}
finishing<-Sys.time()
print(finishing-starting)