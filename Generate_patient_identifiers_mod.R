setwd("/nfs/proj/REPO-TRIAL/montserrat")

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## We generate the identifiers of each patient to generate the patient-specific comorbidity profile ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

if("Identifiers"%in%list.files("Data")==FALSE){dir.create("Data/Identifiers")}
if("Patient_identifiers"%in%list.files("Data/Identifiers")==FALSE){dir.create("Data/Identifiers/Patient_identifiers")}
load("Data/Patient_patient_similarity_network.Rdata")  ;  red000005<-fisher_table

# enfermedad<-unique(gsub("_.+","",unique(c(red000005[,1],red000005[,2]))))

## Remove intra-disease interactions
# uno<-gsub("_.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_.+","",red000005[,2],fixed=F)
uno<-gsub("_GSE.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_GSE.+","",red000005[,2],fixed=F)

tablon<-read.csv2("Patients_information_table_v2.txt",stringsAsFactors = F,sep=",",header=T)
pacientes<-tablon[,1]
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

# Remove diseases already analised from loop:
# Remove from the loop the diseases already analysed
dones <- gsub(".RData","", list.files("Data/Identifiers/Patient_identifiers"))
if(length(dones) > 0) {
  enfermedad <- setdiff(enfermedad,dones)
}



## These are the arguments:
# seq(1,136,28)
# args=commandArgs(trailingOnly = TRUE)  ;  argumento<-as.numeric(args)  ;  print(class(argumento))  ;  if(argumento==113){fin<-23}    if(argumento!=113){fin<-27}
## Start extracting the identifiers
for(b in 1:length(enfermedad)){
  pacientos<-pacientes[grep(enfermedad[b],pacientes, ignore.case = TRUE)]  ;  identificadores<-list()
  for(a in 1:length(pacientos)){
    ini<-Sys.time()
    identificadores[[pacientos[a]]]<-c(grep(pacientos[a],inter000005[,1],ignore.case = TRUE),grep(pacientos[a],inter000005[,2]),ignore.case = TRUE)
    fin<-Sys.time()
    print(paste(a," de ",length(pacientos),sep=""))
  }
  print(paste(enfermedad[b]," finished, ",b," de ",length(enfermedad)," ",fin-ini,sep=""))
  save(identificadores,file=paste("Data/Identifiers/Patient_identifiers/",enfermedad[b],".RData",sep=""))
}
print("Finished!")