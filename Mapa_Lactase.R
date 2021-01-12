tbl4 = read.table(file="//home/victor/LDGH/Diaspora/Output_Admixture/Output_automatiza/K6_Ordenado", header = T, sep = '\t', fill= TRUE)
tbl4_2 <- tbl4[-grep('ID', colnames(tbl4))]
names(tbl4_2)
tbl4_3 <- tbl4_2[-grep('POP', colnames(tbl4_2))]
names(tbl4_3)
nrow <- nrow(tbl4)
nrow



library(rworldmap)

data(countryExData) ### to get the information of countries

LAC=read.table("/home/victor/Projects/Lactose_Intolerance/Mapa/Lactase_frequencies_3.txt", header = TRUE, sep = "\t") #table with ancestry information
#afroDF <- data.frame(country = c("PER", "BRA", "USA","PRI","COL","BRB","CRI","CUB","DOM","VEN","MEX","GTM"),
#                     africa = c(1, 1, 1, 1, 1, 1, 1, 1,1,1,1,1))
afroDF <- data.frame(country = c("PER", "BRA", "PRI","COL","BRB","CRI","CUB","DOM","MEX","GTM"),
                     africa = c(1, 1, 1, 1, 1, 1, 1, 1,1,1))


afroMap <- joinCountryData2Map(afroDF, joinCode = "ISO3",  nameJoinColumn = "country")# This will join your malDF data.frame to the country map data

pdf("//home/victor/Projects/Lactose_Intolerance/Mapa/Mapa_Lactase_ver3.pdf", width = 14.5, height = 12.5)
par(mar = c(1,1,1,1),mgp=c(3,4,4))
#par(fig=c(0,14,3,14)/14)

afroDF <- data.frame(country = c("PER", "BRA", "PRI","COL","BRB","CRI","CUB","DOM","MEX","GTM","USA","URY","CHL"),
                     africa = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))


afroMap <- joinCountryData2Map(afroDF, joinCode = "ISO3",  nameJoinColumn = "country")# This will join your malDF data.frame to the country map data

pdf("//home/victor/Projects/Lactose_Intolerance/Mapa/Mapa_Lactase_Figure3.pdf", width = 14.5, height = 12.5)
par(mar = c(1,1,1,1),mgp=c(3,4,4))
#par(fig=c(0,14,3,14)/14)


#mapParams<-mapCountryData(afroMap, xlim=c(-170,-25), ylim=c(-55,70),nameColumnToPlot="africa", catMethod = "categorical",
mapParams<-mapCountryData(afroMap, xlim=c(-170,-25), ylim=c(-55,70),nameColumnToPlot="africa",numCats=13,
                          missingCountryCol = "aliceblue", addLegend=0, borderCol = "black",oceanCol = "white",colourPalette = c("saddlebrown","olivedrab3","lightsteelblue2",
                                                                                                                                 "royalblue1","orchid4","lightskyblue",
                                                                                                                                 "lightcoral","lightslateblue","firebrick2",
                                                                                                                                 "navy","rosybrown1","orange3",
                                                                                                                                 "yellow1"),
                          aspect=1,mapTitle=" ")

#mapPies(LAC,nameX="LON", nameY="LAT",nameZs=c(names(LAC)[5],names(LAC)[6]),
#                                              zColours=c("palegoldenrod","rosybrown1"),
#                                              xlim=c(-170,-25), ylim=c(-55,70),ratio = 1,
#                                              addCATLegend = "FALSE",symbolSize = 2,maxZVal=1,
#                                              borderCol = "black", oceanCol=NA, landCol=NA,add=TRUE,lty=1)

dev.off()


mapPies(BTW,nameX="LON", nameY="LAT",nameZs=c(names(BTW)[4],names(BTW)[5],names(BTW)[6],names(BTW)[7]),zColours=c("blue","yellow","purple","cyan"), xlim=c(-90,36), ylim=c(-55,50),ratio = 1, addCATLegend = "FALSE",symbolSize = 2.4,maxZVal=1, borderCol = "black", oceanCol=NA, landCol=NA,add=TRUE,lty=1)
arrows(-54, 12, -60, 13, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)  #Barbados
arrows(-47, -40, -54, -30, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Pelotas
arrows(-36, -25, -45, -16, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Bambui
arrows(-32.5, -13, -39.7, -13, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Salvador
arrows(-54, 23, -66, 18, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Puerto Rico
arrows(-121, 29, -116, 36, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#ASW
arrows(-65, 36.5, -85, 36, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#plco
arrows(-83, -13, -76, -13, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Peru
arrows(-93, 5, -76, 5, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Colombia
arrows(53.5, 4, 31.5, 1, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Uganda
arrows(53.5, -10, 38.5, -1, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Kenya
arrows(53.5, -24, 35, -7, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Tanzania
arrows(-21.5, 17.8, -16.5, 13.5, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#Gambia
arrows(-16, 5, -12, 9, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)  #Sierra Leone
arrows(-1, 2, -1, 8, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)  #Ghana
arrows(6, -38.6, 22, -19.2, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#BTW1
arrows(18, -39, 22.5, -19.5, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#BTW2
arrows(34, -38.6, 23, -19.6, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)	#BTW3
arrows(6, 16, 4.5, 7.5, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)  #YRI
arrows(18, 13.5, 6, 8, length = 0, angle = 30, code = 1, col = "black", lty="dashed",lwd=2.5)  #ESN
par(fig=c(0,14,1,3)/14)
par(new=T)
barplot(t(as.matrix(tbl4_3)),yaxt='n',col = c("red1","blue","cyan","yellow","purple","darkgreen"),beside=F,border=NA, xlim=c(0,(nrow+1400)))
legend(x=-560, y = 0.7, c("K=6"), bty = "n",cex=1.5)
par(fig=c(0,14,0,1)/14)
par(new=T)
barplot(t(as.matrix(tbl1_3)),yaxt='n',col = c(0,0, 0,0,0,0,0,0, 0,0,0,0),beside=F,border=NA, xlim=c(0,(nrow+100)))
text(40,0.6,"ESN", las=1,cex=1,srt=90, col="black")
text(160,0.6,"YRI",cex=1,srt=90, col="black")
text(300,0.6,"GWD",cex=1,srt=90, col="black")
text(420,0.6,"MSL",cex=1,srt=90, col="black")
text(900,0.8,"Ghana",cex=1,srt=0, col="black")
text(1780,0.8,"Uganda",cex=1,srt=0, col="black")
text(2080,0.5,"TZWS",cex=1,srt=90, col="black")
text(2160,0.6,"LWK",cex=1,srt=90, col="black")
text(2285,0.23,"BW",cex=1,srt=90, col="black")
text(2290,0.7,"}",cex=2,srt=270, col="black")
text(2253,0.9,"1",cex=0.9,srt=0, col="black")
text(2305,0.9,"2",cex=0.9,srt=0, col="black")
text(2355,0.9,"3",cex=0.9,srt=0, col="black")
text(2430,0.9,"4",cex=0.9,srt=0, col="black")
text(2512,0.9,"5",cex=0.9,srt=0, col="black")
text(2465,0.4,"EUR",cex=0.95,srt=0, col="black")
text(2468,0.7,"}",cex=2.3,srt=270, col="black")
text(2628,0.3,"NAT",cex=0.9,srt=90, col="black")
text(2628,0.75,"}",cex=2,srt=270, col="black")
text(3000,0.6,"PLCO",cex=1,srt=0, col="black")
text(3290,0.6,"ASW",cex=1,srt=90, col="black")
text(3390,0.6,"ACB",cex=1,srt=90, col="black")
text(3460,0.6,"AFP",cex=1,srt=90, col="black")
#arrows(3455,1,3455,0.8,length=0.1,angle=270)
text(3540,0.6,"PUR",cex=1,srt=90, col="black")
text(3605,0.6,"CLM",cex=1,srt=90, col="black")
text(4260,0.8,"Salvador",cex=1,srt=0, col="black")
text(5175,0.8,"Bambui",cex=1,srt=0, col="black")
text(6100,0.8,"Pelotas",cex=1,srt=0, col="black")
dev.off()

