# Patient-specific relative risk estimator
setwd("/nfs/proj/REPO-TRIAL/montserrat/")
# 
load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
red000005<-fisher_table
## Load patient information
tablon<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep=",",header=T)
pacientes<-tablon[,1]
## Remove intra-disease associations
uno<-gsub("_GSE.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_GSE.+","",red000005[,2],fixed=F)
enfermedad<-unique(tablon[,2])

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

intra<-which(uno==dos)  ;  inter000005<-red000005[-intra,]
## Patient identifiers
identidir<-c("Data/Identifiers/Patient_identifiers/")
identifiers<-list()
for(z in 1:length(list.files(identidir))){
  tryCatch({load(paste(identidir,list.files(identidir)[z],sep=""))
  for(q in 1:length(names(identificadores))){identifiers[[names(identificadores)[q]]]<-identificadores[[q]]}
  print(list.files(identidir)[z])}, error = function(e) {print(list.files(identidir[z]))})
}
load("Data/Identifiers/Disease.Rdata") # identificadores
for(a in 1:length(names(identifiers))){
  if(a!=3675){
    print(a)  ;  enf1<-identifiers[[a]]
    for(b in 1:length(enfermedad)){
      if(gsub("_.+","",names(identifiers)[a])!=enfermedad[b]){
        ## Declare objects
        nopresentesiendo<-c() ; presentesiendo<-c() ; presentenosiendo<-c() ; nopresentenosiendo<-c()
        nopresentesiendo1<-c() ; nopresentesiendo2<-c() ; nopresentenosiendo1<-c() ; nopresentenosiendo2<-c()
        ## Select the interactions of interest
        enf2<-identificadores[[enfermedad[b]]] ; par<-intersect(enf1,enf2) ; solo1<-setdiff(enf1,par) ; solo2<-setdiff(enf2,par)
        ## @ @ ##
        ## pRR ##
        siendo<-table(inter000005[par,3])
        # Not present being
        if(length(which(names(siendo)==0))>0){nopresentesiendo1<-as.numeric(siendo[which(names(siendo)==0)])}
        if(length(which(names(siendo)==(-1)))>0){nopresentesiendo2<-as.numeric(siendo[which(names(siendo)==(-1))])}
        if(length(nopresentesiendo1)>0){nopresentesiendo<-nopresentesiendo1}
        if(length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo2}
        if(length(nopresentesiendo1)>0 && length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo1+nopresentesiendo2}
        if(length(nopresentesiendo)==0){nopresentesiendo<-0}                                    ## NEW 22/11/2017
        # Present being
        if(length(which(names(siendo)==1))>0){presentesiendo<-as.numeric(siendo[which(names(siendo)==1)])}
        if(length(presentesiendo)==0){presentesiendo<-0}                                        ## NEW 22/11/2017
        ## Disease
        nosiendo<-table(inter000005[solo1,3])
        # No present not being
        if(length(which(names(nosiendo)==0))>0){nopresentenosiendo1<-as.numeric(nosiendo[which(names(nosiendo)==0)])}
        if(length(which(names(nosiendo)==(-1)))>0){nopresentenosiendo2<-as.numeric(nosiendo[which(names(nosiendo)==(-1))])}
        if(length(nopresentenosiendo1)>0){nopresentenosiendo<-nopresentenosiendo1}
        if(length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo2}
        if(length(nopresentenosiendo1)>0 && length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo1+nopresentenosiendo2}
        if(length(nopresentenosiendo)==0){nopresentenosiendo<-0}                                    ## NEW 22/11/2017
        # Present not being
        if(length(which(names(nosiendo)==1))>0){presentenosiendo<-as.numeric(nosiendo[which(names(nosiendo)==1)])}
        if(length(presentenosiendo)==0){presentenosiendo<-0}                                    ## NEW 22/11/2017
        peventwhenexposed<-presentesiendo/(presentesiendo+nopresentesiendo)
        peventwhennotexposed<-presentenosiendo/(presentenosiendo+nopresentenosiendo)
        division<-peventwhenexposed/peventwhennotexposed
        if(length(division)>0){
          raiz<-sqrt(((nopresentesiendo/presentesiendo)/(nopresentesiendo+presentesiendo))+((nopresentenosiendo/presentenosiendo)/(nopresentenosiendo+presentenosiendo)))
          sup<-exp(log(division)+(1.96*raiz))  ;  inf<-exp(log(division)-(1.96*raiz))
          # RRpos<-rbind(RRpos,c(names(identifiers)[a],enfermedad[b],division,sup,inf))
          jun<-c(names(identifiers)[a],enfermedad[b],division,sup,inf)
          jun<-t(jun)
          write.table(jun,"Data/RR/Patient_pRR.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
        }
        ## @ @ ##
        ## nRR ##
        ## Declare objectes
        nopresentesiendo<-c() ; presentesiendo<-c() ; presentenosiendo<-c() ; nopresentenosiendo<-c()
        nopresentesiendo1<-c() ; nopresentesiendo2<-c() ; nopresentenosiendo1<-c() ; nopresentenosiendo2<-c()
        siendo<-table(inter000005[par,3])
        # Not present being
        if(length(which(names(siendo)==0))>0){nopresentesiendo1<-as.numeric(siendo[which(names(siendo)==0)])}
        if(length(which(names(siendo)==1))>0){nopresentesiendo2<-as.numeric(siendo[which(names(siendo)==1)])}
        if(length(nopresentesiendo1)>0){nopresentesiendo<-nopresentesiendo1}
        if(length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo2}
        if(length(nopresentesiendo1)>0 && length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo1+nopresentesiendo2}
        if(length(nopresentesiendo)==0){nopresentesiendo<-0}                                    ## NEW 22/11/2017
        # Present being
        if(length(which(names(siendo)==(-1)))>0){presentesiendo<-as.numeric(siendo[which(names(siendo)==(-1))])}
        if(length(presentesiendo)==0){presentesiendo<-0}                                    ## NEW 22/11/2017
        ## Disease
        nosiendo<-table(inter000005[solo1,3])
        # Not present not being
        if(length(which(names(nosiendo)==0))>0){nopresentenosiendo1<-as.numeric(nosiendo[which(names(nosiendo)==0)])}
        if(length(which(names(nosiendo)==1))>0){nopresentenosiendo2<-as.numeric(nosiendo[which(names(nosiendo)==1)])}
        if(length(nopresentenosiendo1)>0){nopresentenosiendo<-nopresentenosiendo1}
        if(length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo2}
        if(length(nopresentenosiendo1)>0 && length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo1+nopresentenosiendo2}
        if(length(nopresentenosiendo)==0){nopresentenosiendo<-0}                                    ## NEW 22/11/2017
        # Present not being
        if(length(which(names(nosiendo)==(-1)))>0){presentenosiendo<-as.numeric(nosiendo[which(names(nosiendo)==(-1))])}
        if(length(presentenosiendo)==0){presentenosiendo<-0}                                    ## NEW 22/11/2017
        peventwhenexposed<-presentesiendo/(presentesiendo+nopresentesiendo)
        peventwhennotexposed<-presentenosiendo/(presentenosiendo+nopresentenosiendo)
        division<-peventwhenexposed/peventwhennotexposed
        if(length(division)>0){
          raiz<-sqrt(((nopresentesiendo/presentesiendo)/(nopresentesiendo+presentesiendo))+((nopresentenosiendo/presentenosiendo)/(nopresentenosiendo+presentenosiendo)))
          sup<-exp(log(division)+(1.96*raiz))  ;  inf<-exp(log(division)-(1.96*raiz))
          # RRneg<-rbind(RRneg,c(names(identifiers)[a],enfermedad[b],division,sup,inf))
          jun<-c(names(identifiers)[a],enfermedad[b],division,sup,inf)
          jun<-t(jun)
          write.table(jun,"Data/RR/Patient_nRR.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
        }
      }
    }
  }
}
print("FIN!!")

RRpos<-read.csv2("Data/RR/Patient_pRR.txt",stringsAsFactors = F,sep="\t",header=F) ; RRneg<-read.csv2("Data/RR/Patient_nRR.txt",stringsAsFactors = F,sep="\t",header=F)
RRspositivas<-RRpos[which(as.numeric(RRpos[,5])>1),] ; RRsnegativas<-RRneg[which(as.numeric(RRneg[,5])>1),]
colnames(RRsnegativas)<-c("Patient","Disease","RR","CI_95%_up","CI_95%_down") ; colnames(RRspositivas)<-c("Patient","Disease","RR","CI_95%_up","CI_95%_down")

neg<-RRsnegativas[,1:3] ; pos<-RRspositivas[,1:3]
pasneg<-paste(neg[,1],neg[,2],sep="_") ; paspos<-paste(pos[,1],pos[,2],sep="_")
quit<-intersect(pasneg,paspos)
quitpos<-c();for(a in 1:length(quit)){quitpos<-c(quitpos,which(paspos==quit[a]))}  ;  quitneg<-c();for(a in 1:length(quit)){quitneg<-c(quitneg,which(pasneg==quit[a]))}
negativos<-neg[-quitneg,] ; positivos<-pos[-quitpos,]
po<-cbind(positivos,1) ; colnames(po)[4]<-"Interaction" ; ne<-cbind(negativos,-1) ; colnames(ne)[4]<-"Interaction"
rr1<-rbind(po,ne)
