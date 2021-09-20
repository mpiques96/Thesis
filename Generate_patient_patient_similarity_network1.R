## load the essential packages 
library("affy");library("frma");library("limma");library("genefilter");require(graphics);
library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db");
library("cluster");library("NbClust");library(e1071)

setwd("/nfs/proj/REPO-TRIAL/montserrat")   ## CANVIAR

# Generate required directories
if ("Data"%in%list.files() == FALSE){dir.create("Data")}
if ("ToptableIQR"%in%list.files("Data/") == FALSE){dir.create("Data/ToptableIQR")}
if ("ExpressionSet_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSet_Together")}
if ("ExpressionGenes_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionGenes_Together")}
if ("Correspondences"%in%list.files("Data/") == FALSE){dir.create("Data/Correspondences")}
if ("PatientIQR"%in%list.files("Data/") == FALSE){dir.create("Data/PatientIQR")}
if ("Patienttvalues"%in%list.files("Data/") == FALSE){dir.create("Data/Patienttvalues")}
if ("Networks"%in%list.files() == FALSE){dir.create("Networks")}
if ("Data/ExpressionSetfiles"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSetfiles")}

#Added by me:
if ("Data/Barcodes"%in%list.files("Data/") == FALSE){dir.create("Data/Barcodes")}


## set seed ##
## @@ @@ @@ ##
set.seed(1)

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
##                       Generate gene expression matrices for each disease separated by studies                      ##
## conduct differential expression analyses comparing in each study all the case samples with all the control samples ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

rawdatadir<-"/nfs/scratch/m.selles/raw_data/"   ### MOURE DATA
rawdatafiles<-list.files(rawdatadir)
rawdatafiles

## load gene symbol - probe association and the probes that must be removed 
load("Data/Remove_Symbols_u133plus2.Rdata") # remo
load("Data/Gene_Symbols_u133plus2.Rdata") # simbolos

## indicate which diseases have been already analyzed 
if("Analyzed_diseases.txt"%in%list.files("Data/") == TRUE){echos<-read.table("Data/Analyzed_diseases.txt",stringsAsFactors = F)[,1];rawdatafiles<-setdiff(rawdatafiles,echos);dones<-echos}
if("Analyzed_diseases.txt"%in%list.files("Data/") == FALSE){dones<-c()}
print("We start the loop!")

for (a in 1:length(rawdatafiles)){
  inicio<-Sys.time()
  ## print information regarding the dataset that is going to be analyzed
  print(paste(rawdatafiles[a],"Control:",length(list.files(paste(rawdatadir,rawdatafiles[a],"/Control/",sep=""))),"Case:",
              length(list.files(paste(rawdatadir,rawdatafiles[a],"/Case/",sep=""))),sep="      "))
  control<-list.files(paste(rawdatadir,"/",rawdatafiles[a],"/Control",sep=""))
  path_con<-paste(rawdatadir,rawdatafiles[a],"/Control/",control,sep="")
  case<-list.files(paste(rawdatadir,"/",rawdatafiles[a],"/Case",sep=""))
  path_ca<-paste(rawdatadir,rawdatafiles[a],"/Case/",case,sep="")
  path_read<-c(path_con, path_ca)
  controls<-cbind(control,rep("control",length(control)))
  cases<-cbind(case,rep("case",length(case)))
  targets<-rbind(controls, cases)
  affyBatch<-ReadAffy(filenames=path_read)  ###PROBLEMA
  expSet<-frma(affyBatch)
  save(expSet,file=paste("Data/ExpressionSetfiles/",rawdatafiles[a],".Rdata",sep=""))
  barcoded<-barcode(expSet)
  save(barcoded,file=paste("Data/Barcodes/",rawdatafiles[a],".Rdata",sep=""))
  comoda<-t(as.data.frame(expSet))
  ## save the normalized expression matrix (with probes)
  save(comoda,file=paste("Data/ExpressionSet_Together/",rawdatafiles[a],".Rdata",sep=""))
  comoda1<-comoda[-length(comoda[,1]),]
  ## remove probes not associated to gene symbols
  comoda2<-comoda1[-remo,]
  coluna<-colnames(comoda2)
  genes<-unique(simbolos)
  expression<-c()
  ## convert the probes' expression matrix into a gene expression matrix, median values are calculated when seveeral probes refer to the same gene symbol
  for (c in 1:length(genes)){
    com<-which(simbolos==as.character(genes[c]))
    if (length(com)>1){
      iq<-apply(comoda[com,],2,median)
    }
    if (length(com)==1){
      iq<-comoda[com,]
    }
    expression<-rbind(expression,iq)
  }
  colnames(expression)<-coluna
  rownames(expression)<-genes
  columna<-c(control,case)
  ## change samples names
  newpatientnames<-c(paste("Control_",1:length(control),sep=""),paste("Patient_",1:length(case),sep=""))
  correspondencia<-cbind(columna,newpatientnames)
  colnames(correspondencia)<-c("Original_name","Patient_name")
  ## save patient name - sample ID information
  save(correspondencia,file=paste("Data/Correspondences/",rawdatafiles[a],".Rdata",sep=""))
  colnames(expression)<-newpatientnames
  controls<-cbind(paste("Control_",1:length(control),sep=""),rep("control",length(control)))
  cases<-cbind(paste("Patient_",1:length(case),sep=""),rep("case",length(case)))
  targets<-rbind(controls, cases)
  ## save gene expression matrix with the new names
  save(expression,file=paste("Data/ExpressionGenes_Together/",rawdatafiles[a],".Rdata",sep=""))
  ## conduct differential gene expression analysis using LIMMA
  design<-cbind(CONTROL=c(rep(1,length(control)),rep(0,length(case))), CASE=c(rep(0,length(control)),rep(1,length(case))))
  rownames(design)<-targets[,1]
  cont.matrix<-makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
  ## get DEGs from IQR
  fit<-lmFit(expression,design)
  fit2<-contrasts.fit(fit, cont.matrix)
  fit2<-eBayes(fit2)
  toptableIQR<-topTable(fit2, number=length(fit$Amean), adjust.method="BH", sort.by="p")
  ## save all patient vs. all control differential expression results
  write.table(toptableIQR,paste("Data/ToptableIQR/",rawdatafiles[a],".txt",sep=""),sep="\t",quote=F)
  print(paste(a,"  out of  ",length(rawdatafiles),":  ",rawdatafiles[a],sep=""))
  dones<-c(dones,rawdatafiles[a])
  ## save the name of the analyzed disease (with its' corresponding study)
  write.table(dones,"Data/Analyzed_diseases.txt",quote=F,col.names=F,row.names=F,sep="\t")
  final<-Sys.time()
  diferencia<-final-inicio
  print(diferencia)
}
print("Gene-expression matrices and global differential expression analyses finished!!")