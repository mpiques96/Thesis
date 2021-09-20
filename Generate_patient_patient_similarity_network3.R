## load the essential packages 
library("affy");library("frma");library("limma");library("genefilter");require(graphics);
library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db");
library("cluster");library("NbClust");library(e1071)

setwd("/nfs/proj/REPO-TRIAL/montserrat")


## Extract top X up- and down-regulated genes and generate a matrix of tvalues for LINCS analysis ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## duration: 35.88323 mins
topini<-Sys.time()
directory<-"Data/Patienttvalues/"
## indicate the number of DEG we want to select
seleccion<-500
patients<-list.files(directory)
nombres<-gsub(".Rdata","",patients,fixed=T)
## avoid re-analyzing patients
if (paste("Top_",seleccion,"_sDEGs.Rdata",sep="")%in%list.files("Data/") == TRUE){
  load(paste("Data/Top_",seleccion,"_sDEGs.Rdata",sep=""))
  hechos<-paste(names(sdegs),".Rdata",sep="")
  echos<-names(sdegs)
  patients<-setdiff(patients,hechos)
  nombres<-setdiff(nombres,echos) ; print("We do have sDEGs")
}
if (paste("Top_",seleccion,"_sDEGs.Rdata",sep="")%in%list.files("Data/") == FALSE){sdegs<-list() ; print("We do not have sDEGs")}
## select the up- and down-regulated genes
if(length(patients)>0){
  load(paste(directory,patients[1],sep=""))
  genes<-names(tvalues)
  tablatvalores<-c()
  for(a in 1:length(patients)){
    # a<-1
    load(paste(directory,patients[a],sep=""))
    tablatvalores<-cbind(tablatvalores,as.numeric(tvalues[genes]))
    tvalores<-sort(tvalues)
    sdegs[[nombres[a]]]$down<-tvalores[1:seleccion]
    sdegs[[nombres[a]]]$up<-tvalores[(length(tvalores)-(seleccion-1)):length(tvalores)]
    print(paste(a,"  out of  ",length(patients),"  patients",sep=""))
  }
  rownames(tablatvalores)<-genes
  colnames(tablatvalores)<-nombres
  save(sdegs,file=paste("Data/Top_",seleccion,"_sDEGs.Rdata",sep=""))
  saveRDS(tablatvalores,file="Data/LINCS_tvalues_patients.rds")
  print(paste("Top ",seleccion,"genes selected!!",sep="  "))
}
topfin<-Sys.time()
print(topfin-topini)



## Generate a list with DEGs names ##
## @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ ##
## duration: 1.329571 mins
indeg<-Sys.time()
patients<-list.files(directory)
nombres<-gsub(".Rdata","",patients,fixed=T)
## avoid re-analyzing patients
if (paste("Top_",seleccion,"_sDEGs_names.Rdata",sep="")%in%list.files("Data/") == FALSE){degs<-list() ; print("We do not have the list with sDEGs names")}
if (paste("Top_",seleccion,"_sDEGs_names.Rdata",sep="")%in%list.files("Data/") == TRUE){
  load(paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))
  hechos<-paste(names(degs),".Rdata",sep="")
  echos<-names(degs)
  patients<-setdiff(patients,hechos)
  nombres<-setdiff(nombres,echos)
  print("We have the list with sDEGs names")
}
## select the names of the up- and down-regulated genes
if(length(patients)>0){
  for(a in 1:length(patients)){
    load(paste(directory,patients[a],sep=""))
    tvalores<-sort(tvalues)
    degs[[nombres[a]]]$down<-names(tvalores[1:seleccion])
    degs[[nombres[a]]]$up<-names(tvalores[(length(tvalores)-(seleccion-1)):length(tvalores)])
    print(paste(a,"  out of  ",length(patients),"  finished",sep=""))
  }
  save(degs,file=paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))
}
fideg<-Sys.time()
print(fideg-indeg)
fisin<-Sys.time()
load(paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))


## Fisher test analysis ##
## @@ @@ @@ @@ @@ @@ @@ ##
## duration: 17.61 hours
background<-20502   # wtf?
load("Data/Top_500_sDEGs_names.Rdata")
fisin<-Sys.time()
todos<-c()
toditos<-c()
for (a in 1:(length(names(degs))-1)){
  positi1<-degs[[a]]$up ; negati1<-degs[[a]]$down ; pega<-c() ; pegados<-c()
  for (b in (a+1):length(names(degs))){
    positi2<-degs[[b]]$up ; negati2<-degs[[b]]$down
    interspd <- length(intersect(positi1,positi2))
    kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
    fispd<-fisher.test(kkpd,alternative="greater") $p.value
    intersnd <- length(intersect(negati1,negati2))
    kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
    fisnd<-fisher.test(kknd,alternative="greater") $p.value
    interspi <- length(intersect(positi1,negati2))
    kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
    fispi<-fisher.test(kkpi,alternative="greater") $p.value
    intersni <- length(intersect(negati1,positi2))
    kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
    fisni<-fisher.test(kkni,alternative="greater") $p.value
    junt<-c(names(degs)[a],names(degs)[b],fispi,fispd,fisni,fisnd)
    jun<-t(junt)
    write.table(jun,"Fisher_overlaps.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
  }
  print(a/(length(names(degs))-1))
}
print("Finished!!") ; fisfin<-Sys.time() ; print(fisfin-fisin)

## read the file and generate the network
fis<-read.csv("Fisher_overlaps.txt",stringsAsFactors = F,sep="\t",header=F)
threshold<-0.001
negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
## save the patient-patient similarity network including 0s
fisher_table<-cbind(fis[,1:2],interactions) ; colnames(fisher_table)<-c("Disease_1","Disease_2","Interaction")
save(fisher_table,file="Data/Patient_patient_similarity_network.Rdata")
network<-fisher_table[which(fisher_table[,3]!=0),]
## save the generated "patient-patient similarity network"
write.table(network,"Networks/Patient_patient_similarity_network.txt",sep="\t",quote=F,row.names=F)

## Generate a matrix of 1s, 0s and -1s with patients in rows and genes in columns ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 2.74 hours
starting<-Sys.time()
patients<-list.files("Data/Patienttvalues/") ; nome<-gsub(".Rdata","",patients,fixed=T) ; onesandzerosm<-list()
for(a in 1:length(patients)){
  load(paste("Data/Patienttvalues/",patients[a],sep=""))
  updown<-rep(0,length(tvalues)) ; updown[1:500]<-(-1) ; updown[(length(updown)-(500-1)):length(updown)]<-1 ; names(updown)<-names(tvalues) ; onesandzerosm[[nome[a]]]<-updown
}
genes<-names(onesandzerosm[[1]])
newmatrix<-c()
for(a in 1:length(names(onesandzerosm))){newmatrix<-rbind(newmatrix,onesandzerosm[[a]][genes]) ; print(paste(a," out of ",length(names(onesandzerosm)),sep=""))}
rownames(newmatrix)<-nome
save(newmatrix,file="Data/Matrix_with_patients_and_genes_1_0_-1.Rdata")
finishing<-Sys.time() ; print(finishing-starting)

print("End of part 1")