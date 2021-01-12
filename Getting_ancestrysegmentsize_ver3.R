

### ACB and ASW
msp.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/1000Genomes/ACB_ASW_continental_chr2.msp.tsv"

msp.file<-read.table(msp.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)


##look for the interval
spos_snp<-findInterval(135851076,msp.file$spos) ### to look for the start position SNP of the region

#ncol(msp.file)

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(msp.file))){                             ## scannin thought all individuals
  if (msp.file[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(msp.file[i])

    while (msp.file[spos_snp-beginP,i]==2) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (msp.file[spos_snp+endP,i]==2) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    startpositions [inds.with.target.ancestry]<- msp.file[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- msp.file[spos_snp+endP,2]
  }
}

pop.asw.acb.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.asw.acb.size$V1<-ids
pop.asw.acb.size$V2<-startpositions
pop.asw.acb.size$V3<-endpositions

for (i in c(1:length(ids))){    
  pop.asw.acb.size$V4[i]<-round((pop.asw.acb.size$V3[i]-pop.asw.acb.size$V2[i])/1000000,3)
}
library(dplyr)


pop.asw.acb.size$V5<-gsub("HG[0-9]+\\.[0-9]", 'ACB', pop.asw.acb.size$V1)
pop.asw.acb.size$V5<-gsub("NA[0-9]+\\.[0-9]", 'ASW', pop.asw.acb.size$V5)


names(pop.asw.acb.size)[names(pop.asw.acb.size) == "V1"] <- "ID"
names(pop.asw.acb.size)[names(pop.asw.acb.size) == "V2"] <- "spos"
names(pop.asw.acb.size)[names(pop.asw.acb.size) == "V3"] <- "epos"
names(pop.asw.acb.size)[names(pop.asw.acb.size) == "V4"] <- "size"
names(pop.asw.acb.size)[names(pop.asw.acb.size) == "V5"] <- "Pop"

#ACB<-pop.LA.size %>%filter( pop.LA.size$V1 %in% grep("NA",pop.LA.size$V1,value=TRUE))
#mean(ACB$V4)
#sd(ACB$V4)

#Peruvian Genome Project

msp.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/PGP/PGP_continental_chr2.msp.tsv"
pgp.LA<-read.table(msp.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

##look for the interval
spos_snp<-findInterval(135851076,pgp.LA$spos) ### to look for the start position SNP of the region

library(dplyr)
#pgp.ids<-subset(LPids.file %>%filter( LPids.file$V1 %in% grep("pgp",LPids.file$V1,value=TRUE)),select = 2)
#pgp.position<-grep(paste(pgp.ids$V2, collapse="|"),colnames(LPLA.file))
#pgp.LA<-subset(LPLA.file,select = c(1:6,pgp.position ))
spos_snp<-findInterval(135851076,pgp.LA$spos) ### to look for the start position SNP of the region

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(pgp.LA))){                             ## scannin thought all individuals
  if (pgp.LA[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(pgp.LA[i])

    while (isTRUE(pgp.LA[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(pgp.LA[spos_snp+endP,i]==2) && spos_snp+endP<nrow(pgp.LA)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    
    startpositions [inds.with.target.ancestry]<- pgp.LA[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- pgp.LA[spos_snp+endP,2]
  }
}

pop.pgp.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.pgp.size$V1<-ids
pop.pgp.size$V2<-startpositions
pop.pgp.size$V3<-endpositions

for (i in c(1:length(ids))){    
  pop.pgp.size$V4[i]<-round((pop.pgp.size$V3[i]-pop.pgp.size$V2[i])/1000000,3)
}

pop.pgp.size$V5<-gsub('[0-9]+.[0-9]', '', pop.pgp.size$V1)   ## To remove number from the column in order to create population names

Ayacucho<-pop.pgp.size %>%filter( pop.pgp.size$V1 %in% grep("Ayacucho",pop.pgp.size$V1,value=TRUE))
mean(Ayacucho$V4)
sd(Ayacucho$V4)


pop.pgp.size2<-pop.pgp.size[!(pop.pgp.size$V5=="Aimaras" | pop.pgp.size$V5=="Chachapoyas" | pop.pgp.size$V5=="Chopccas" | pop.pgp.size$V5=="Jacarus"
                              | pop.pgp.size$V5=="Moche" | pop.pgp.size$V5=="Quechuas" | pop.pgp.size$V5=="Shipibo" 
                              | pop.pgp.size$V5=="Tallanes" | pop.pgp.size$V5=="Uros"),]


names(pop.pgp.size2)[names(pop.pgp.size2) == "V1"] <- "ID"
names(pop.pgp.size2)[names(pop.pgp.size2) == "V2"] <- "spos"
names(pop.pgp.size2)[names(pop.pgp.size2) == "V3"] <- "epos"
names(pop.pgp.size2)[names(pop.pgp.size2) == "V4"] <- "size"
names(pop.pgp.size2)[names(pop.pgp.size2) == "V5"] <- "Pop"


### LARGE PD

largepdids.path<-"/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/LARGE-PD/ids_threecolumns.txt"
largepdLA.path<-"/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/LARGE-PD/large-PD.local_ancestry.3_groups.chr2.msp.tsv"
LPids.file<-read.table(largepdids.path, header = F, sep = '', fill= TRUE)
LPLA.file<-read.table(largepdLA.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

spos_snp<-findInterval(135851076,LPLA.file$spos) ### to look for the start position SNP of the region

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(LPLA.file))){                             ## scannin thought all individuals
  if (LPLA.file[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(LPLA.file[i])

    while (isTRUE(LPLA.file[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(LPLA.file[spos_snp+endP,i]==2) && spos_snp+endP<nrow(LPLA.file)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    
    startpositions [inds.with.target.ancestry]<- LPLA.file[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- LPLA.file[spos_snp+endP,2]
  }
}

pop.LA.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.LA.size$V1<-ids
pop.LA.size$V2<-startpositions
pop.LA.size$V3<-endpositions

for (i in c(1:length(ids))){    
  pop.LA.size$V4[i]<-round((pop.LA.size$V3[i]-pop.LA.size$V2[i])/1000000,3)
}

pop.LA.size$V5<-gsub("\\.[0-9]", '', pop.LA.size$V1)

pop.large.ids<-subset(merge(pop.LA.size, LPids.file, by.x = "V5", by.y = "V3", all.x = TRUE), select = c(1,3,4,5,6))

names(pop.large.ids)[names(pop.large.ids) == "V5"] <- "ID"
names(pop.large.ids)[names(pop.large.ids) == "V2.x"] <- "spos"
names(pop.large.ids)[names(pop.large.ids) == "V3"] <- "epos"
names(pop.large.ids)[names(pop.large.ids) == "V4"] <- "size"
names(pop.large.ids)[names(pop.large.ids) == "V1.y"] <- "Pop"

pop.large.ids2<-pop.large.ids[!(pop.large.ids$Pop=="Peru" | pop.large.ids$Pop=="Peru_Puno" | pop.large.ids$Pop=="Colombia_Medellin") ,]

library("tidyverse")
library("ggplot2")

large_pgp<-rbind(pop.pgp.size2,pop.large.ids2)

#### EPIGEN

msp.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/EPIGEN/epigen_continental_chr2.msp.tsv"
epigen.LA<-read.table(msp.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

##look for the interval

library(dplyr)
spos_snp<-findInterval(135851076,epigen.LA$spos) ### to look for the start position SNP of the region

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(epigen.LA))){                             ## scannin thought all individuals
  if (epigen.LA[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(epigen.LA[i])

    while (isTRUE(epigen.LA[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(epigen.LA[spos_snp+endP,i]==2) && spos_snp+endP<nrow(epigen.LA)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    
    startpositions [inds.with.target.ancestry]<- epigen.LA[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- epigen.LA[spos_snp+endP,2]
  }
}

pop.epigen.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.epigen.size$V1<-ids
pop.epigen.size$V2<-startpositions
pop.epigen.size$V3<-endpositions

Bambui.haps<-pop.epigen.size[grep("Bambui", pop.epigen.size$V1), ]
removebambuihap<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Bambui_list_haps2remove", header = FALSE)
Bambui.norelated<-filter(Bambui.haps, !V1%in%removebambuihap$V1) ## remove
Pelotas.haps<-pop.epigen.size[grep("Pelotas", pop.epigen.size$V1), ]
removepelotashap<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Pelotas_list_haps2remove", header = FALSE)
Pelotas.norelated<-filter(Pelotas.haps, !V1%in%removepelotashap$V1) ## remove

Salvador.haps<-pop.epigen.size[grep("SCAALA", pop.epigen.size$V1), ]
removesalvadorhap<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Salvador_list_haps2remove", header = FALSE)
Salvador.norelated<-filter(Salvador.haps, !V1%in%removesalvadorhap$V1) ## remove




for (i in c(1:length(ids))){    
    pop.epigen.size$V4[i]<-round((pop.epigen.size$V3[i]-pop.epigen.size$V2[i])/1000000,3)
}

pop.epigen.size$V5<-gsub('_[a-zA-Z][0-9]+\\.+[0-9]', '', pop.epigen.size$V1)
pop.epigen.size$V5<-gsub('_B[0-9]+r\\.+[0-9]', '', pop.epigen.size$V5)
pop.epigen.size$V5<-gsub('_S[0-9]+B\\.+[0-9]', '', pop.epigen.size$V5)

names(pop.epigen.size)[names(pop.epigen.size) == "V1"] <- "ID"
names(pop.epigen.size)[names(pop.epigen.size) == "V2"] <- "spos"
names(pop.epigen.size)[names(pop.epigen.size) == "V3"] <- "epos"
names(pop.epigen.size)[names(pop.epigen.size) == "V4"] <- "size"
names(pop.epigen.size)[names(pop.epigen.size) == "V5"] <- "Pop"

## Guatemala

msp.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/Guatemala/Guatemala_continental_chr2.msp.tsv"
Guatemala.LA<-read.table(msp.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

##look for the interval

library(dplyr)
spos_snp<-findInterval(135851076,Guatemala.LA$spos) ### to look for the start position SNP of the region

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(Guatemala.LA))){                             ## scannin thought all individuals
  if (Guatemala.LA[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(Guatemala.LA[i])
    
    while (isTRUE(Guatemala.LA[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(Guatemala.LA[spos_snp+endP,i]==2) && spos_snp+endP<nrow(Guatemala.LA)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    
    startpositions [inds.with.target.ancestry]<- Guatemala.LA[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- Guatemala.LA[spos_snp+endP,2]
  }
}

pop.Guatemala.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.Guatemala.size$V1<-ids
pop.Guatemala.size$V2<-startpositions
pop.Guatemala.size$V3<-endpositions

for (i in c(1:length(ids))){    
    pop.Guatemala.size$V4[i]<-round((pop.Guatemala.size$V3[i]-pop.Guatemala.size$V2[i])/1000000,3)
}

pop.Guatemala.size$V5<-pop.Guatemala.size[pop.Guatemala.size$V5]<-"Guatemala"

names(pop.Guatemala.size)[names(pop.Guatemala.size) == "V1"] <- "ID"
names(pop.Guatemala.size)[names(pop.Guatemala.size) == "V2"] <- "spos"
names(pop.Guatemala.size)[names(pop.Guatemala.size) == "V3"] <- "epos"
names(pop.Guatemala.size)[names(pop.Guatemala.size) == "V4"] <- "size"
names(pop.Guatemala.size)[names(pop.Guatemala.size) == "V5"] <- "Pop"


#### MXL CLM and PUR

msp.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/1000Genomes/MXL_CLM_PUR_continental_chr2.msp.tsv"
msp.file<-read.table(msp.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

##look for the interval
spos_snp<-findInterval(135851076,msp.file$spos) ### to look for the start position SNP of the region

#ncol(msp.file)

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(msp.file))){                             ## scannin thought all individuals
  if (msp.file[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(msp.file[i])
    
    ### backward
    while (isTRUE(msp.file[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(msp.file[spos_snp+endP,i]==2)&& spos_snp+endP<nrow(msp.file)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    startpositions [inds.with.target.ancestry]<- msp.file[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- msp.file[spos_snp+endP,2]
  }
}

pop.3pop.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.3pop.size$V1<-ids
pop.3pop.size$V2<-startpositions
pop.3pop.size$V3<-endpositions

for (i in c(1:length(ids))){    
  
  pop.3pop.size$V4[i]<-round((pop.3pop.size$V3[i]-pop.3pop.size$V2[i])/1000000,3)
  
}
library(dplyr)

pop.3pop.size$V5<-gsub("_[a-zA-Z]+[0-9]+\\.+[0-9]", '', pop.3pop.size$V1)

names(pop.3pop.size)[names(pop.3pop.size) == "V1"] <- "ID"
names(pop.3pop.size)[names(pop.3pop.size) == "V2"] <- "spos"
names(pop.3pop.size)[names(pop.3pop.size) == "V3"] <- "epos"
names(pop.3pop.size)[names(pop.3pop.size) == "V4"] <- "size"
names(pop.3pop.size)[names(pop.3pop.size) == "V5"] <- "Pop"

#### TOPMED

topmed.path<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/TOPMED/Central_America_continental_chr2.msp.tsv"
topmed.file<-read.table(topmed.path, header = T, sep = '\t', fill= TRUE,comment.char = '&', skip=1)

spos_snp<-findInterval(135851076,topmed.file$spos) ### to look for the start position SNP of the region

inds.with.target.ancestry<-0
ids<-c()
startpositions <- c()
endpositions <- c()

for (i in c(7:ncol(topmed.file))){                             ## scannin thought all individuals
  if (topmed.file[spos_snp,i]==2){                             ## identifying individuals with EUR ancestry in the interval region
    inds.with.target.ancestry<-inds.with.target.ancestry+1
    beginP<-0
    endP<-0
    ids [inds.with.target.ancestry]<- colnames(topmed.file[i])
    
    ### backward
    while (isTRUE(topmed.file[spos_snp-beginP,i]==2) && spos_snp-beginP>1) {                ## identifying the start position of the region
      beginP<-beginP+1
    }
    while (isTRUE(topmed.file[spos_snp+endP,i]==2)&& spos_snp+endP<nrow(topmed.file)) {                  ## identifying the end postion of the region
      endP<-endP+1
    }
    startpositions [inds.with.target.ancestry]<- topmed.file[spos_snp-beginP,3]
    endpositions [inds.with.target.ancestry]<- topmed.file[spos_snp+endP,2]
  }
}

pop.top.size <- as.data.frame(matrix(ncol=4,nrow=length(startpositions)))

pop.top.size$V1<-ids
pop.top.size$V2<-startpositions
pop.top.size$V3<-endpositions

for (i in c(1:length(ids))){    
    pop.top.size$V4[i]<-round((pop.top.size$V3[i]-pop.top.size$V2[i])/1000000,3)
}

library(dplyr)

pop.top.size$V5<-gsub("_[0-9]+\\.+[0-9]", '', pop.top.size$V1)

names(pop.top.size)[names(pop.top.size) == "V1"] <- "ID"
names(pop.top.size)[names(pop.top.size) == "V2"] <- "spos"
names(pop.top.size)[names(pop.top.size) == "V3"] <- "epos"
names(pop.top.size)[names(pop.top.size) == "V4"] <- "size"
names(pop.top.size)[names(pop.top.size) == "V5"] <- "Pop"






#### Merging

large_pgp_epigen.1kgp<-rbind(large_pgp,pop.epigen.size,pop.asw.acb.size,pop.Guatemala.size,pop.3pop.size,pop.top.size)

large_pgp_epigen.1kgp$Pop<- factor(large_pgp_epigen.1kgp$Pop , levels=c("ASW","ACB","MXL","Guatemala","CostaRican", "Cuban","Dominican","PUR",
                                                                        "Colombia_Bogota","CLM","Afro_des", "Tumbes", "Lambayeque", "Ancash","Trujillo","Lima", 
                                                                        "Arequipa", "Moquegua", "Tacna","Ayacucho", "Cusco", "Puno", "Iquitos","Chile", 
                                                                        "SCAALA", "Bambui", "Brazil_SaoPaolo", "Brazil_RibeiraoPreto", "Brazil_PortoAlegre", "Pelotas",
                                                                        "Uruguay"),ordered = TRUE)


colores<-c("rosybrown1", "lightgrey", "firebrick2","navy","lightskyblue","yellow1","lightslateblue","lightsteelblue2",
           "royalblue1","royalblue1", "deeppink","saddlebrown","saddlebrown","saddlebrown","saddlebrown","saddlebrown",
           "saddlebrown","saddlebrown","saddlebrown","saddlebrown","saddlebrown","saddlebrown","saddlebrown","magenta3",
           "olivedrab3","olivedrab3","olivedrab3", "olivedrab3","olivedrab3","olivedrab3",
           "orange3")
#ggplot(large_pgp_epigen.1kgp, aes(x=Pop, y=size)) + geom_boxplot(fill=colores, color="black")



#write.table(large_pgp_epigen.1kgp,"/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/Getting_Hap_size_results.txt",
#            col.names=T,row.names=F,quote =F)


library("tidyverse")
library("ggplot2")
library(dplyr)
large_pgp_epigen.1kgp<-read.table("/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/Getting_Hap_size_results.txt",header = TRUE)
epigen2remove<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/EPIGEN2remove",header = FALSE)

Cuban.haps<-large_pgp_epigen.1kgp[grep("Cuba", large_pgp_epigen.1kgp$Pop), ]
Dominican.haps<-large_pgp_epigen.1kgp[grep("Dom", large_pgp_epigen.1kgp$Pop), ]
keepcubahap<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Cuban3.txt", header = FALSE)
Cuba.remove<-filter(Cuban.haps, !ID%in%keepcubahap$V1) ## remove
keepdomhap<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Dominican3.txt", header = FALSE)
Dom.remove<-filter(Dominican.haps, !ID%in%keepdomhap$V1) ## remove

removeinds.dom<-as.character(Dom.remove$ID)
removeinds.cuba<-as.character(Cuba.remove$ID)
removeinds.bras<-as.character(epigen2remove$V1)


individuals2remove<-c(removeinds.bras,removeinds.cuba,removeinds.dom)
#keepcubahap<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Cuban3.txt", header = FALSE)
#Cuba.norelated<-filter(Cuban.haps, !ID%in%keepcubahap$V1) ## remove
#keepdomhap<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Dominican3.txt", header = FALSE)
#Dom.norelated<-filter(Dominican.haps, !ID%in%keepdomhap$V1) ## remove

large_pgp_epigen.1kgp.unrelated<-filter(large_pgp_epigen.1kgp, !ID%in%individuals2remove) ## remove

large_pgp_epigen.1kgp.unrelated$Pop<- factor(large_pgp_epigen.1kgp.unrelated$Pop , levels=c("ASW","ACB","MXL","Guatemala","CostaRican", "Cuban","Dominican","PUR",
                                                                        "Colombia_Bogota","CLM","Afro_des", "Tumbes", "Lambayeque", "Ancash","Trujillo","Lima", 
                                                                        "Arequipa", "Moquegua", "Tacna","Ayacucho", "Cusco", "Puno", "Iquitos","Chile", 
                                                                        "SCAALA", "Bambui", "Brazil_SaoPaolo", "Brazil_RibeiraoPreto", "Brazil_PortoAlegre", "Pelotas",
                                                                        "Uruguay"),ordered = TRUE)




graph_latin<-ggplot(large_pgp_epigen.1kgp.unrelated, aes(x=Pop, y=size)) + geom_boxplot(fill=colores, color="black")

graph_latin<-  graph_latin + theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5),axis.text.y = element_text(angle = 90, size = 8, vjust = 0.5))


graph_latin<-graph_latin+ scale_x_discrete(name= "", breaks=c("ASW","ACB","MXL","Guatemala","CostaRican", "Cuban","Dominican","PUR",
                                       "Colombia_Bogota","CLM","Afro_des", "Tumbes", "Lambayeque", "Ancash","Trujillo","Lima", 
                                       "Arequipa", "Moquegua", "Tacna","Ayacucho", "Cusco", "Puno", "Iquitos","Chile", 
                                       "SCAALA", "Bambui", "Brazil_SaoPaolo", "Brazil_RibeiraoPreto", "Brazil_PortoAlegre", "Pelotas",
                                       "Uruguay"),
                              labels=c("Afro-Americans", "Afro-Caribbeans", "Mexico", "Guatemala", "Costa Rica", "Cuba", "Dominican Republic", "Puerto Rico",
                                       "Colombia (Bogota)","Colombia (Medellin)","Afro-Peruvians","Peru (Tumbes)","Peru (Lambayeque)","Peru (Ancash)","Peru (Trujillo)","Peru (Lima)",
                                       "Peru (Arequipa)","Peru (Moquegua)","Peru (Tacna)","Peru (Ayacucho)","Peru (Cusco)","Peru (Puno)","Peru (Iquitos)","Chile",
                                       "Brazil (Salvador)","Brazil (Bambui)","Brazil (Sao Paulo)","Brazil (Ribeirao Preto)","Brazil (Porto Alegre)","Brazil (Pelotas)",
                                       "Uruguay"))
#graph_latin<-graph_latin + scale_y_discrete(name= "size (Mb)",breaks=c("0","50","100","150","200","250"),limits=c("0","50","100","150","200","250"))
graph_latin
ggsave("/home/victor/Projects/Lactose_Intolerance/Local_ancestry/Results/December_boxplot_size_ver2021.jpg",scale = 1.5, width = 14, height = 7, units = "cm",dpi = 600)


dev.off()
mean(pop.LA.size$V4)
sd(pop.LA.size$V4)