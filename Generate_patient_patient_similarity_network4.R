## load the essential packages 
library("affy");library("frma");library("limma");library("genefilter");require(graphics);
library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db");
library("cluster");library("NbClust");library(e1071)

setwd("/nfs/proj/REPO-TRIAL/montserrat")


## @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ ##
## Calculate the number of intra-disease interactions increasing the number of selected genes ##
## @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ ##

if("Supplementary_documents"%in%list.files()==FALSE){dir.create("Supplementary_documents")}
if("Different_number_of_genes"%in%list.files("./Supplementary_documents")==FALSE){dir.create("Supplementary_documents/Different_number_of_genes")}
if("Plots"%in%list.files()==FALSE){dir.create("Plots")}

## Run lines from 334 to 387 assigning the following numbers to the object "selec" --> 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000
starting<-Sys.time()
patients<-list.files("Data/Patienttvalues/") ; nome<-gsub(".Rdata","",patients,fixed=T) ; degs<-list()
## we select the number of genes as differentially expressed (100,200,300,400,500,1000,2000,3000,4000,5000)
####### MODIFIED #########
number2 <- c(100,200,300,400,500,1000,2000)
for(number in number2){
  for(a in 1:length(patients)){
    load(paste("Data/Patienttvalues/",patients[a],sep=""))
    tvalues<-tvalues[order(as.numeric(tvalues),decreasing=T)]
    up<-names(tvalues)[1:number] ; down<-names(tvalues)[(length(tvalues)-(number-1)):length(tvalues)]
    degs[[gsub(".Rdata","",patients[a])]]$up<-up ; degs[[gsub(".Rdata","",patients[a])]]$down<-down
  }
  ## Fisher test ##
  background<-20502
  diseases<-gsub("_.+","",names(degs)) ; disease<-unique(diseases)
  for(z in 1:length(disease)){
    cojo<-which(diseases==disease[z])
    digs<-list() ; for(y in cojo){digs[[names(degs)[y]]]<-degs[y]}
    for (a in 1:(length(names(digs))-1)){
      positi1<-digs[[a]][[1]]$up ; negati1<-digs[[a]][[1]]$down ; pega<-c() ; pegados<-c()
      for (b in (a+1):length(names(digs))){
        positi2<-digs[[b]][[1]]$up ; negati2<-digs[[b]][[1]]$down
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
        junt<-c(names(digs)[a],names(digs)[b],fispi,fispd,fisni,fisnd)
        jun<-t(junt)
        write.table(jun,paste("Supplementary_documents/Different_number_of_genes/Fisher_overlaps_",number,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F,append = T)
      }
      print(a/(length(names(digs))-1))
    }
  }
  ficheros<-list.files("Supplementary_documents/Different_number_of_genes/")
  suma<-c()
  for(a in 1:length(ficheros)){
    ## read the file and generate the network
    fis<-read.csv(paste("Supplementary_documents/Different_number_of_genes/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t",header=F) ; threshold<-0.05
    negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
    poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
    interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
    tabint<-table(interactions) ; neg<-as.numeric(tabint[which(names(tabint)=="-1")]) ; pos<-as.numeric(tabint[which(names(tabint)=="1")]) ; suma<-c(suma,neg+pos)
  }
  seleccion<-as.numeric(gsub(".txt","",gsub("Fisher_overlaps_","",ficheros)))
  pdf(file="Plots/Intra-disease_interactions_varying_number_of_deg.pdf")
  plot(seleccion,suma,xlab="Number of Differentially Expressed Genes selected",ylab="Number of intra-disease interactions")
  dev.off()
}

###### END OF MODIFICATION #######