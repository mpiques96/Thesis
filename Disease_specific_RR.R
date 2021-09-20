setwd("/nfs/proj/REPO-TRIAL/montserrat/")
## Extract identifiers and then calculate Relative Risk between ICD9s, ICD10s, diseases & patient-subgroups ##
## Confidence intervals ##
## 95% --> 1.96
## 99% --> 2.576
## load function
estimate_rr<-function(disease,identificadores,inter0001,path,confidenceinterval){
  RRpos<-c()
  RRneg<-c()
  for(a in 1:(length(disease)-1)){
    beginloop<-Sys.time()
    ## calculate relative risks for each disease 
    print(disease[a])
    ## select first diseases' interactions of interest
    enf1<-identificadores[[disease[a]]]
    for(b in (a+1):length(disease)){
      notpresentbeing<-c() ; presentbeing<-c() ; presentnotbeing<-c() ; notpresentnotbeing<-c() ; notpresentbeing1<-c() ; notpresentbeing2<-c()
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c()
      ## select second diseases' interactions of interest
      enf2<-identificadores[[disease[b]]]
      par<-intersect(enf1,enf2)
      only1<-setdiff(enf1,par)
      only2<-setdiff(enf2,par)
      
      ## @ @ ##
      ## pRR ##
      ## @ @ ##
      being<-table(inter0001[par,3])
      # Not present being
      if(length(which(names(being)==0))>0){notpresentbeing1<-as.numeric(being[which(names(being)==0)])}
      if(length(which(names(being)==(-1)))>0){notpresentbeing2<-as.numeric(being[which(names(being)==(-1))])}
      if(length(notpresentbeing1)>0){notpresentbeing<-notpresentbeing1}
      if(length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing2}
      if(length(notpresentbeing1)>0 && length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing1+notpresentbeing2}
      if(length(notpresentbeing)==0){notpresentbeing<-0}
      # Present being
      if(length(which(names(being)==1))>0){presentbeing<-as.numeric(being[which(names(being)==1)])}
      if(length(presentbeing)==0){presentbeing<-0}
      ## Disease 1 ##
      ## @ @@ @@ @ ##
      notbeing<-table(inter0001[only1,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==(-1)))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==1))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRpos<-rbind(RRpos,c(disease[a],disease[b],division,sup,inf))
        # print(paste(disease[a],disease[b],division,sup,inf,sep="  "))
      }
      ## Disease 2 ##
      ## @ @@ @@ @ ##
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c() ; presentnotbeing<-c() ; notbeing<-table(inter0001[only2,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==(-1)))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==1))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRpos<-rbind(RRpos,c(disease[b],disease[a],division,sup,inf))
        # print(paste(disease[b],disease[a],division,sup,inf,sep="  "))
      }
      
      ## @ @ ##
      ## nRR ##
      ## @ @ ##
      notpresentbeing<-c() ; presentbeing<-c() ; presentnotbeing<-c() ; notpresentnotbeing<-c() ; notpresentbeing1<-c() ; notpresentbeing2<-c()
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c()
      being<-table(inter0001[par,3])
      # Not present being
      if(length(which(names(being)==0))>0){notpresentbeing1<-as.numeric(being[which(names(being)==0)])}
      if(length(which(names(being)==1))>0){notpresentbeing2<-as.numeric(being[which(names(being)==1)])}
      if(length(notpresentbeing1)>0){notpresentbeing<-notpresentbeing1}
      if(length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing2}
      if(length(notpresentbeing1)>0 && length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing1+notpresentbeing2}
      if(length(notpresentbeing)==0){notpresentbeing<-0}
      # Present being
      if(length(which(names(being)==(-1)))>0){presentbeing<-as.numeric(being[which(names(being)==(-1))])}
      if(length(presentbeing)==0){presentbeing<-0}
      ## Disease 1 ##
      ## @ @@ @@ @ ##
      notbeing<-table(inter0001[only1,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==1))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==(-1)))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRneg<-rbind(RRneg,c(disease[a],disease[b],division,sup,inf))
        # print(paste(disease[a],disease[b],division,sup,inf,sep="  "))
      }
      ## Disease 2 ##
      ## @ @@ @@ @ ##
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c() ; presentnotbeing<-c()
      notbeing<-table(inter0001[only2,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==1))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==(-1)))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRneg<-rbind(RRneg,c(disease[b],disease[a],division,sup,inf))
        # print(paste(disease[b],disease[a],division,sup,inf,sep="  "))
      }
    }
    endloop<-Sys.time()
    save(RRpos,file=paste(path,"pRR.Rdata",sep="_"))
    save(RRneg,file=paste(path,"nRR.Rdata",sep="_"))
    print(paste("Finished ",a," out of ",length(disease)," in ",endloop-beginloop,sep=" "))
  }
  rrs<-list("pRR"=RRpos,"nRR"=RRneg)
  print("Finished!!")
  return(rrs)
}

## generate folders
if("RR"%in%list.files("Data/")==FALSE){dir.create("./Data/RR/")}
if ("Identifiers"%in%list.files("Data/") == FALSE){dir.create("Data/Identifiers")}



## @@ @@ @@ ##
## Diseases ##
## @@ @@ @@ ##
## Generate identifiers ##
## @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 34.12 mins
beggining<-Sys.time()
load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
red000005<-fisher_table
## Load patient information
tablon<-read.csv2("Patients_information_table_v2.txt",stringsAsFactors = F,sep=",",header=T)
## Remove intra-disease associations
uno<-gsub("_GSE.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_GSE.+","",red000005[,2],fixed=F)
disease<-unique(tablon[,2])

map1 <- c("Asthma_central","Asthma_peripheral","Atopic_dermatitis_non_lesional","Endometriosis_early_mild","Endometriosis_early_severe",
          "Endometriosis_mid_mild","Endometriosis_mid_severe","Endometriosis_prol_mild","Endometriosis_prol_severe","Ependymoma_1",
          "Ependymoma_2","Erythematotelangiectatic_rosacea","Human_Papillomavirus_Cervix_1","Human_Papillomavirus_Cervix_2",
          "Obesity_fat","Obesity_skin","Papulopustular_rosacea","Phymatous_rosacea","Stage_1_Huntingtons_disease","Stage_2_Huntingtons_disease",
          "Vitiligo_lesional","Vitiligo_non_lesional","Vitiligo_peri_lesional")
map2 <- c("Asthma","Asthma","Atopic_dermatitis","Endometriosis","Endometriosis","Endometriosis","Endometriosis","Endometriosis",
          "Endometriosis", "Ependymoma","Ependymoma","Rosacea","Human_Papillomavirus_Cervix","Human_Papillomavirus_Cervix",
          "Obesity","Obesity","Rosacea","Rosacea","Huntingtons_disease","Huntingtons_disease","Vitiligo","Vitiligo","Vitiligo")
mapping <- data.frame(map1,map2)

for(i in 1:length(map1)){
  uno[uno %in% mapping$map1[i]] = mapping$map2[i]
  dos[dos %in% mapping$map1[i]] = mapping$map2[i]
}

intra<-which(uno==dos)  ;  inter0001<-red000005[-intra,]  ;  identificadores<-list()
print("Starting the loop")
for(a in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[a]]]<-c(grep(paste(disease[a],"_",sep=""),inter0001[,1]),grep(paste(disease[a],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(a," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file="Data/Identifiers/Disease.Rdata")
ending<-Sys.time()
print(paste("Finished in: ",ending-beggining,"!",sep=""))
## Calculate Relative Risks ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
## estimated duration: 12.35 minutes
starting<-Sys.time()
## load previously generated identifiers
load("Data/Identifiers/Disease.Rdata")
path<-"./Data/RR/Disease"
estimate_rr(disease,identificadores,inter0001,path,1.96)
finishing<-Sys.time()  ;  print(finishing-starting)