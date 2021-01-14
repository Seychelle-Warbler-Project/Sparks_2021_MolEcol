# # -----------------------------------------------------------------------

# ----- Telomere heritability manuscript
# ----- AMS
# ----- PAC/MAC models
# ----- Power analysis 

# # -----------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(lme4)
library(simr)
library(ggplot2)


#~~ data 

telodata <- read.csv("Sparks_Telo_h2_MS_Data_2020.csv", header=T)


#~~ set factors etc

colnames(telodata)[which(names(telodata) == "BirdID")] <- "ID"

telodata$ID <- as.factor(telodata$ID)
telodata$Sex <- as.factor(telodata$Sex)
telodata$FPID <- as.factor(telodata$FPID)
telodata$mum <- as.factor(telodata$mum)
telodata$dad <- as.factor(telodata$dad)
telodata$U_PlateID <- as.factor(telodata$U_PlateID)
telodata$Technician <- as.factor(telodata$Technician)
telodata$sqrtRTL <- sqrt(telodata$RTL)
telodata$logAge <- log10(telodata$AgeY)


# # -----------------------------------------------------------------------

#~~ Power analysis

# # -----------------------------------------------------------------------


#~~ make model subset

telodata.all<-droplevels(subset(telodata, !is.na(AgeY) & !is.na(Sex) & !is.na(Technician) & !is.na(MAC) & !is.na(PAC)))


#~~ z-transform sqrt RTL values

telodata.all$cRTL <- (telodata.all$sqrtRTL - mean(telodata.all$sqrtRTL, na.rm=TRUE)) / sd(telodata.all$sqrtRTL, na.rm=TRUE)


#~~ PAC/MAC model

PACmodel <- lmer(cRTL ~ logAge + Sex + Technician + MAC + PAC +
                    (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                  data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACmodel)


#~~ run power analysis - part 1 - effect sizes 0.001-0.025 

results <- data.frame(successes=NA, trials=NA, mean=NA, lower=NA, upper=NA)

for(i in seq(0.001, 0.025, 0.001)) {
  fixef(PACmodel)['PAC']<- i
  temp<-summary(powerSim(PACmodel, fixed("PAC", "z"), nsim=500))
  results <- rbind(results, temp)
}


results$PACeff <- c(NA, seq(0.001, 0.025, 0.001))
results <- results[-1,]

write.csv(results, file="PACpowerresults190220_1.csv")


#~~ run power analysis - part 2 - effect sizes 0.026-0.050 

rm(results)

results <- data.frame(successes=NA, trials=NA, mean=NA, lower=NA, upper=NA)

for(i in seq(0.026, 0.050, 0.001)) {
  fixef(PACmodel)['PAC']<- i
  temp<-summary(powerSim(PACmodel, fixed("PAC", "z"), nsim=500))
  results <- rbind(results, temp)
}


results$PACeff <- c(NA, seq(0.026, 0.050, 0.001))
results <- results[-1,]

write.csv(results, file="PACpowerresults190220_2.csv")


# # -----------------------------------------------------------------------

#~~ PAC power analysis graph for supplementary material - figure s5

# # -----------------------------------------------------------------------


#~~ call in both results files and merge

results1 <- read.csv("PACpowerresults190220_1.csv", header=T)
results2 <- read.csv("PACpowerresults190220_2.csv", header=T)

PowerData <- rbind(results1, results2)

#~~ converting regression slope (beta) estimate from lme4 to a correlation coefficient 
#~~ using r = beta slope*(sd x / sd y)

0.02*(sd(telodata.all$PAC)/sd(telodata.all$cRTL))

#~~ graph code

PACpowerplot <- ggplot(PowerData, aes(x=PACeff,y=mean)) +
  labs(x="PAC effect size", y="Power") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=10, colour="grey26")) +
  geom_pointrange(aes(ymin=lower, ymax=upper), width=.2) +
  geom_segment(aes(x = 0.02, y = 0, xend = 0.02, yend = 0.8), linetype="dashed", color = "grey66") +
  geom_segment(aes(x = 0, y = 0.8, xend = 0.02, yend = 0.8), linetype="dashed", color = "grey66") + 
  geom_point(shape=19, colour="grey12", size=1)

PACpowerplot

tiff("FigureS5_PACPowerGraph.tif", width=2500, height=2000, res=400)

PACpowerplot

dev.off()


