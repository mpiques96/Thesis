## load the essential packages 
library("affy");library("frma");library("limma");library("genefilter");require(graphics);
library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db");
library("cluster");library("NbClust");library(e1071)

setwd("/nfs/proj/REPO-TRIAL/montserrat")

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Conduct differential expression analyses comparing each case sample with all the control samples from the same study ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

expresiondir<-c("Data/ExpressionGenes_Together/")
expresionfiles<-list.files(expresiondir)
nombrecitos<-gsub(".Rdata","",expresionfiles,fixed = T)
## avoid re-analyzing patients
if(length(list.files("Data/Patienttvalues/"))>0){
  pacientes<-list.files("Data/Patienttvalues/")
  if(length(pacientes)>0){
    hecho<-c()
    dones<-gsub("_Patient.+",".Rdata",pacientes)
    done<-unique(dones)
    for(a in done){
      load(paste(expresiondir,a,sep=""))
      if(length(which(dones==a))==length(grep("Patient",colnames(expression)))){
        hecho<-c(hecho,a)
      }
    }
  }
  if(length(pacientes)==0){
    hecho<-c()
  }
  expresionfiles<-setdiff(expresionfiles,hecho)
}
## start patient-specific differential expression analysis
for(a in 1:length(expresionfiles)){
  load(paste(expresiondir,expresionfiles[a],sep=""))
  patientes<-grep("Patient",colnames(expression))
  controles<-grep("Control",colnames(expression))
  for(b in 1:length(patientes)){
    expresion<-expression[,c(controles,patientes[b])]
    control<-colnames(expression)[controles]
    case<-colnames(expression)[patientes[b]]
    controls<-cbind(paste("Control_",1:length(control),sep=""),rep("control",length(control)))
    cases<-cbind(paste("Patient_",1:length(case),sep=""),rep("case",length(case)))
    targets<-rbind(controls, cases)
    design<-cbind(CONTROL=c(rep(1,length(control)),rep(0,length(case))), CASE=c(rep(0,length(control)),rep(1,length(case))))
    rownames(design)<-targets[,1]
    cont.matrix<-makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
    fit<-lmFit(expresion,design)  ##getting DEGs from IQR
    fit2<-contrasts.fit(fit, cont.matrix)
    fit2<-eBayes(fit2)
    toptableIQR<-topTable(fit2, number=length(fit$Amean), adjust.method="BH", sort.by="p")
    write.table(toptableIQR,paste("./Data/PatientIQR/",nombrecitos[a],"_",case,".txt",sep=""),sep="\t",quote=F)
    tvalues<-toptableIQR[,3]
    names(tvalues)<-rownames(toptableIQR)
    save(tvalues,file=paste("./Data/Patienttvalues/",nombrecitos[a],"_",case,".Rdata",sep=""))
  }
  print(paste(expresionfiles[a],"  we have finished  ",a,"  out of  ",length(expresionfiles),sep=""))
}
print("Patient-specific differential gene expression finished!!")