

Peru.path<- "/home/victor/Projects/Lactose_Intolerance/Peru_frequencies.txt"
Peru.freqs<-read.table(Peru.path, header = T, sep = '\t', fill= TRUE)

library(ggplot2)

# Init Ggplot
g<-ggplot(Peru.freqs, aes(x=European_proportion,y=Allelic_percentage)) + geom_point()
g<-g + xlim(c(0, 20)) + ylim(c(0, 30))
g + geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)

hap <-ggplot(Peru.freqs, aes(x=European_proportion,y=haplotype_percentage)) + geom_point()
hap<-hap + xlim(c(0, 30)) + ylim(c(0, 30))
hap + geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)

#===================================================================================

wwf.path<- "/home/victor/Projects/Lactose_Intolerance/WorldWide_frequencies.txt"

wwf.freqs<-read.table(wwf.path, header = T, sep = '\t', fill= TRUE)

library(ggplot2)
library(gridExtra)

# Init Ggplot
g<-ggplot(wwf.freqs, aes(x=Mean_European_proportion,y=Allelic_percentage)) +
  geom_point(aes(colour=as.character(Country)),size=3.5) +
  scale_colour_manual(name="Populations",
                      breaks=as.character(wwf.freqs$Country),
                      values=as.character(wwf.freqs$Color)) +
    xlim(c(0, 100)) + ylim(c(0, 100))+labs(x="Average of European ancestry (%)", y="-13910*T frequency (%)") +
    geom_abline(intercept = 0, slope = 1, color="firebrick2", linetype="dashed",  size=0.5, alpha=0.4) +
    theme(text = element_text(size = 16))

jpeg(filename = "/home/victor/Projects/Lactose_Intolerance/December_Allelic_frequency_vs_EUR_ver2.jpg",width = 20, height = 15, units = "cm", pointsize =12,  res = 600)
g
dev.off()

hap <-ggplot(wwf.freqs, aes(x=Mean_European_proportion,y=haplotype_percentage)) + 
    geom_point(aes(colour=as.character(Country)),size=3.5) +
    scale_colour_manual(name="Populations",
                      breaks=as.character(wwf.freqs$Country),
                      values=as.character(wwf.freqs$Color)) +
    xlim(c(0, 100)) + ylim(c(0, 100))+labs(x="Average of European ancestry (%)", y="European Haplotype frequency (%)") +
    geom_abline(intercept = 0, slope = 1, color="firebrick2", linetype="dashed", size=0.5, alpha=0.4) +
    theme(text = element_text(size = 16))

jpeg(filename = "/home/victor/Projects/Lactose_Intolerance/December_HAP_frequency_vs_EUR_ver2.jpg",width = 20, height = 15, units = "cm", pointsize =12,  res = 600)
hap
dev.off()


##### EPIGEN ADMIXTURE
library(dplyr)

epigen.file<- "/home/victor/Projects/Lactose_Intolerance/admixture_results/EPIGEN_admixture.proportions"
epigen.admixture<-read.table(epigen.file, header = FALSE, sep = "", fill= TRUE) ## "" to any number of space in the table

colnames(epigen.admixture) <- c("ID","NAT", "EUR","AFR","EAS")

Bambui.info<-epigen.admixture[grep("Bambui", epigen.admixture$ID), ]
removebambui<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Bambui_list", header = FALSE)
Bambui.info.norelated<-filter(Bambui.info, !ID%in%removebambui$V1)
mean(Bambui.info.norelated$EUR)

Pelotas.info<-epigen.admixture[grep("Pelotas", epigen.admixture$ID), ]
removepelotas<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Pelotas_list", header = FALSE)
Pelotas.info.norelated<-filter(Pelotas.info, !ID%in%removepelotas$V1)
mean(Pelotas.info.norelated$EUR)

Salvador.info<-epigen.admixture[grep("SCAALA", epigen.admixture$ID), ]
removesalvador<-read.table("/home/victor/Projects/Lactose_Intolerance/EPIGEN/Salvador_list", header = FALSE)
Salvador.info.norelated<-filter(Salvador.info, !ID%in%removesalvador$V1)
mean(Salvador.info.norelated$EUR)


topmed.file<- "/home/victor/Projects/Lactose_Intolerance/admixture_results/3pops_topmed_admixture.proportions"
topmed.admixture<-read.table(topmed.file, header = FALSE, sep = "", fill= TRUE) ## "" to any number of space in the table

colnames(topmed.admixture) <- c("ID","NAT", "EUR","AFR","EAS")

CR.info<-topmed.admixture[grep("Costa", topmed.admixture$ID), ]
mean(CR.info$EUR)

Cuba.info<-topmed.admixture[grep("Cuba", topmed.admixture$ID), ]
keepcuba<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Cuban2.txt", header = FALSE)
Cuba.info.norelated<-filter(Cuba.info, ID%in%keepcuba$V1) ## KEEP
mean(Cuba.info.norelated$EUR)

Dom.info<-topmed.admixture[grep("Dom", topmed.admixture$ID), ]
keepdom<-read.table("/home/victor/Projects/Lactose_Intolerance/TOPMED/IDs_HCHS_SOL_Dominican2.txt", header = FALSE)
Dom.info.norelated<-filter(Dom.info, ID%in%keepdom$V1) ## KEEP
mean(Dom.info.norelated$EUR)


###Local ancestry

large.file<- "/home/victor/Projects/Lactose_Intolerance/Local_ancestry/three_columns_ID_HAP_ancestry.txt"
large.la<-read.table(large.file, header = FALSE, sep = "", fill= TRUE) ## "" to any number of space in the table

large.id<-read.table("/home/victor/Projects/Lactose_Intolerance/Local_ancestry/largepd_old_newfam_id.txt", header = FALSE, sep = "", fill= TRUE) ## "" to any number of space in the table

full.large<-merge(large.id, large.la, by.x = "V4", by.y = "V1", all.y = T)
small.large<-subset(full.large,select = c(4:6))

Uruguay.info<-small.large[grep("Uruguay", small.large$V3.x), ]

head(Uruguay.info)


nchar(gsub("European", "", Uruguay.info$V3.y))

lengths(lapply(Uruguay.info$V3.y, grepRaw, pattern = "European", all = TRUE, fixed = TRUE))

as.data.frame(table(Uruguay.info$V3.y))


Chile.info<-small.large[grep("Chile", small.large$V3.x), ]
as.data.frame(table(Chile.info$V3.y))


Ribeirao.info<-small.large[grep("Ribeirao", small.large$V3.x), ]
as.data.frame(table(Ribeirao.info$V3.y))

SP.info<-small.large[grep("Paolo", small.large$V3.x), ]
as.data.frame(table(SP.info$V3.y))


##### CORRELATION

tables4<- "/home/victor/Projects/Lactose_Intolerance/Supplementary_table_4_alt.txt"
tableS4<-read.table(tables4, header = T, sep = '\t', fill= TRUE)

library(dplyr)
tableS4_EURmore<-tableS4[(tableS4$Mean_European_proportion>50),]
tableS4_EURless<-tableS4[(tableS4$Mean_European_proportion<50),]


res1<-cor.test(tableS4$rs4988235.T,tableS4$Mean_European_proportion, method = "pearson")
res1

res1<-cor.test(tableS4$EUR_haplotype_frequency,tableS4$Mean_European_proportion, method = "pearson")
res1



res1<-cor.test(tableS4$rs4988235.T,tableS4$Mean_European_proportion, method = "pearson")
res1

reslower<-cor.test(tableS4_EURless$rs4988235.T,tableS4_EURless$Mean_European_proportion, method = "pearson")
reslower

resmore<-cor.test(tableS4_EURmore$rs4988235.T,tableS4_EURmore$Mean_European_proportion, method = "pearson")
resmore


library(ggplot2)

#wwf.path<- "/home/victor/Projects/Lactose_Intolerance/WorldWide_frequencies.txt"

#wwf.freqs<-read.table(wwf.path, header = T, sep = '\t', fill= TRUE)

library(ggplot2)
library(gridExtra)

# Init Ggplot
g<-ggplot(tableS4, aes(x=Mean_European_proportion,y=rs4988235.T*100)) +
  geom_point(aes(colour=as.character(Country)),size=3.5) +
  scale_colour_manual(name="Populations",
                      breaks=as.character(tableS4$Country),
                      values=as.character(tableS4$Color)) +
  xlim(c(0, 100)) + ylim(c(0, 100))+labs(x="Average of European ancestry (%)", y="-13910*T frequency (%)") +
  geom_abline(intercept = 0, slope = 1, color="firebrick2", linetype="dashed",  size=0.5, alpha=0.4) +
  theme(text = element_text(size = 16))

jpeg(filename = "/home/victor/Projects/Lactose_Intolerance/December_Allelic_frequency_vs_EUR_2021.jpg",width = 20, height = 15, units = "cm", pointsize =12,  res = 600)
g
dev.off()

hap <-ggplot(tableS4, aes(x=Mean_European_proportion,y=EUR_haplotype_frequency*100)) + 
  geom_point(aes(colour=as.character(Country)),size=3.5) +
  scale_colour_manual(name="Populations",
                      breaks=as.character(tableS4$Country),
                      values=as.character(tableS4$Color)) +
  xlim(c(0, 100)) + ylim(c(0, 100))+labs(x="Average of European ancestry (%)", y="European Haplotype frequency (%)") +
  geom_abline(intercept = 0, slope = 1, color="firebrick2", linetype="dashed", size=0.5, alpha=0.4) +
  theme(text = element_text(size = 16))

jpeg(filename = "/home/victor/Projects/Lactose_Intolerance/December_HAP_frequency_vs_EUR_ver3.jpg",width = 20, height = 15, units = "cm", pointsize =12,  res = 600)
hap
dev.off()


g<-ggplot(tableS4, aes(x=EUR_haplotype_frequency*100,y=proportion.Allele.EURprop*100)) +
  geom_point(aes(colour=as.character(Country)),size=3.5) +
  scale_colour_manual(name="Populations",
                      breaks=as.character(tableS4$Country),
                      values=as.character(tableS4$Color)) +
  xlim(c(0, 100)) + ylim(c(0, 100))+labs(x="European Haplotype frequency (%)", y="Proportion of -13910*T in European haplotypes (%)") +
  geom_abline(intercept = 0, slope = 1, color="firebrick2", linetype="dashed",  size=0.5, alpha=0.4) +
  theme(text = element_text(size = 16))

jpeg(filename = "/home/victor/Projects/Lactose_Intolerance/December_Allelic_proportion_vs_EUR_Haplotypes.jpg",width = 20, height = 15, units = "cm", pointsize =12,  res = 600)
g
dev.off()


