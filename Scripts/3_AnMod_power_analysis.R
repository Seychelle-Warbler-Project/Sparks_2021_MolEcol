# # -----------------------------------------------------------------------

# ----- Telomere heritability manuscript
# ----- AMS
# ----- Animal models
# ----- Power analysis

# # -----------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(pedantics)
library(MasterBayes)
source("phensim.R") # from Michael Morrissey 
library(ggplot2)


#~~ read in pruned ped for RTL analysis

ped <- read.csv("Telo_h2_Pruned_Pedigree_2020.csv", header=T)


#~~ order ped - list mums/dads before offspring

pedtelo_order <- orderPed(ped) 


#~~ telomere h2 power analysis

n<-1000
heritabilities<-seq(0,0.5,by=0.01)
p_values<-matrix(NA,n,length(heritabilities))
power<-array(dim=length(heritabilities))

for(x in 1:length(heritabilities)) {
  for(y in 1:n) {
    simPhen<-phensim(pedtelo_order,randomA=heritabilities[x], randomE=1-heritabilities[x])$phenotypes
    simPhen$mum<-pedtelo_order$mum[match(simPhen$id,pedtelo_order$id)]
    simPhen$mumPhenotype<-simPhen$trait_1[match(simPhen$mum,as.numeric(as.character(simPhen$id)))]
    p_values[y,x]<-anova(lm(simPhen$trait_1~simPhen$mumPhenotype))[1,5]}
  power[x]<-table(p_values[,x]<0.05)["TRUE"]/n 
}


#~~ save power results

power_results<-cbind(heritabilities,power)
power_results
write.table(power_results, file="Teloh2_Ped_h2_power_250220.csv", row.names = F, sep = ",", quote = F)


#~~ graph for supplementary material - power analysis graph - figure s7

PowerData <- read.csv("Teloh2_Ped_h2_power_250220.csv", header=T)

pedpowerplot <- ggplot(PowerData, aes(x=heritabilities,y=power)) +
  labs(x="Heritability", y="Power") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=10, colour="grey26")) +
  geom_segment(aes(x = 0.17, y = 0, xend = 0.17, yend = 0.805), linetype="dashed", color = "grey66") +
  geom_segment(aes(x = 0, y = 0.805, xend = 0.17, yend = 0.805), linetype="dashed", color = "grey66") + 
  geom_point(shape=19, colour="grey12", size=1)

pedpowerplot

tiff("FigureS7_PedPowerGraph.tif", width=2500, height=2000, res=400)

pedpowerplot

dev.off()


# # -----------------------------------------------------------------------

# Power analysis for the discussion

# # -----------------------------------------------------------------------


#~~ power analysis using chick data only

#~~ read in telomere data

telodata <- read.csv("Sparks_Telo_h2_MS_Data_2020.csv", header=T)

colnames(telodata)[which(names(telodata) == "BirdID")] <- "ID"

telodatach <- droplevels(subset(telodata, AgeClass=="CH"))
telochicks <- telodatach$ID

keepIDs2<-unique(telochicks) # 319 birds


#~~ prune pedigree keeping only informative individuals for this analysis

pedteloch<-prunePed(ped, keep=keepIDs2) # retains 802 birds


#~~ order ped 

pedteloch_order <- orderPed(pedteloch) # list mums/dads before offspring


#~~ telomere h2 power analysis

n<-1000
heritabilities<-seq(0,0.5,by=0.01)
p_values<-matrix(NA,n,length(heritabilities))
power<-array(dim=length(heritabilities))

for(x in 1:length(heritabilities)) {
  for(y in 1:n) {
    simPhen<-phensim(pedteloch_order,randomA=heritabilities[x], randomE=1-heritabilities[x])$phenotypes
    simPhen$mum<-pedteloch_order$mum[match(simPhen$id,pedteloch_order$id)]
    simPhen$mumPhenotype<-simPhen$trait_1[match(simPhen$mum,as.numeric(as.character(simPhen$id)))]
    p_values[y,x]<-anova(lm(simPhen$trait_1~simPhen$mumPhenotype))[1,5]}
  power[x]<-table(p_values[,x]<0.05)["TRUE"]/n 
}


#~~ save power results

power_results_ch<-cbind(heritabilities,power)
power_results_ch
write.table(power_results_ch, file="Teloh2_Ped_h2_power_ch_180420.csv", row.names = F, sep = ",", quote = F)


#~~ power for birds >7 years

telodatasen <- droplevels(subset(telodata, AgeY>7))
telosen <- telodatasen$ID

keepIDs3<-unique(telosen) # 161 birds


#~~ prune pedigree keeping only informative individuals for this analysis

pedtelosen<-prunePed(ped, keep=keepIDs3) # retains 364 birds


#~~ order ped 

pedtelosen_order <- orderPed(pedtelosen) # list dams/sires before offspring


#~~ telomere h2 power analysis

n<-1000
heritabilities<-seq(0,0.5,by=0.01)
p_values<-matrix(NA,n,length(heritabilities))
power<-array(dim=length(heritabilities))

for(x in 1:length(heritabilities)) {
  for(y in 1:n) {
    simPhen<-phensim(pedtelosen_order,randomA=heritabilities[x], randomE=1-heritabilities[x])$phenotypes
    simPhen$mum<-pedtelosen_order$mum[match(simPhen$id,pedtelosen_order$id)]
    simPhen$mumPhenotype<-simPhen$trait_1[match(simPhen$mum,as.numeric(as.character(simPhen$id)))]
    p_values[y,x]<-anova(lm(simPhen$trait_1~simPhen$mumPhenotype))[1,5]}
  power[x]<-table(p_values[,x]<0.05)["TRUE"]/n 
}


#~~ save power results

power_results<-cbind(heritabilities,power)
power_results
write.table(power_results, file="Teloh2_Ped_h2_power_sen_180420.csv", row.names = F, sep = ",", quote = F)
