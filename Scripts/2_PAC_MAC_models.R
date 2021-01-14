# # -----------------------------------------------------------------------

# ----- Telomere heritability manuscript
# ----- AMS
# ----- PAC/MAC models
# ----- The models 

# # -----------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(ggplot2)
library(lme4)
library(dplyr)
library(cowplot)
library(data.table)
library(AICcmodavg)


#~~ data 

telodata <- read.csv("Sparks_Telo_h2_MS_Data_2020.csv", header=T)


#~~ function - variance inflation factors


vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


# # -----------------------------------------------------------------------

# MAC/PAC summary stats and graphs for supplementary material 

# # -----------------------------------------------------------------------


#~~ make unique offspring dataset for MAC/PAC stats

telodataoffspring <- telodata[ ! duplicated (telodata$BirdID), ]
telodataoffspring$mum <- as.factor(telodataoffspring$mum)
telodataoffspring$dad <- as.factor(telodataoffspring$dad)
telodataoffspring$BirdID <- as.factor(telodataoffspring$BirdID)
telodataoffspring <- droplevels(telodataoffspring)

#~~ correlation

cor.test(telodataoffspring$MAC, telodataoffspring$PAC)

#~~ summary stats

summary(telodataoffspring$MAC)
summary(telodataoffspring$PAC)

#~~ number of mums and dads with multiple offspring

table(table(telodataoffspring$dad))
table(table(telodataoffspring$mum))


#~~ graph for supplementary material - figure s6 - histograms of MAC/PAC and correlation between MAC and PAC

MAChist<-ggplot(telodataoffspring, aes(x=MAC)) +
  geom_histogram(alpha=.5, position="identity") +
  xlab("Maternal age at conception") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey26")) 
PAChist<-ggplot(telodataoffspring, aes(x=PAC)) +
  geom_histogram(alpha=.5, position="identity") +
  xlab("Paternal age at conception") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey26")) 
MACPACcor<-ggplot(telodataoffspring,aes(x=MAC,y=PAC)) + geom_point(shape=21, alpha=0.5, colour="grey4")+
  labs(x="Maternal age at conception", y="Paternal age at conception") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey26")) 

plot_grid(MAChist, PAChist, MACPACcor, labels = "AUTO", ncol = 3, nrow=1)

plot<-plot_grid(MAChist, PAChist, MACPACcor, labels = "AUTO", ncol = 3, nrow=1)

tiff("FigureS6_PACMACGraphs.tif", width=4000, height=1250, res=400)

plot

dev.off()


# # -----------------------------------------------------------------------

# MAC/PAC models

# # -----------------------------------------------------------------------


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

# PAC/MAC across all ages

# # -----------------------------------------------------------------------


#~~ model subset

telodata.all<-droplevels(subset(telodata, !is.na(AgeY) & !is.na(Sex) & !is.na(Technician) & !is.na(MAC) & !is.na(PAC)))


#~~ z-transform sqrt RTL values

telodata.all$cRTL <- (telodata.all$sqrtRTL - mean(telodata.all$sqrtRTL, na.rm=TRUE)) / sd(telodata.all$sqrtRTL, na.rm=TRUE)


#~~ make mean MAC/PAC and deviation from mean MAC/PAC

telodata.all<-telodata.all %>% group_by(mum) %>% mutate(MeanMAC = mean(MAC)) %>% ungroup()
telodata.all$DevMeanMAC<-telodata.all$MAC-telodata.all$MeanMAC

telodata.all<-telodata.all %>% group_by(dad) %>% mutate(MeanPAC = mean(PAC)) %>% ungroup()
telodata.all$DevMeanPAC<-telodata.all$PAC-telodata.all$MeanPAC


#~~ model

PACall <- lmer(cRTL ~ logAge + Sex + Technician + MAC + PAC +
                    (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                  data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACall)


#~~ check diagnostics

plot(PACall)
qqnorm(resid(PACall))
vif.mer(PACall)


#~~ save model output

PACallmodel<-as.data.frame((summary(PACall)$coefficients))
PACallmodel<-setDT(PACallmodel, keep.rownames = TRUE)[]

write.csv(PACallmodel, file="PACMACall190220.csv")

PACallmodel2<-as.data.frame(VarCorr(PACall))

write.csv(PACallmodel2, file="PACMACallranef190220.csv")


#~~ check sig of fixed effects

#~~ age - x=94.767      df=1  p<2.2e-16 

PACall1<-update(PACall, ~ . -logAge)
anova(PACall1, PACall)


#~~ sex - x=0.31      df=1     p=0.5777

PACall2<-update(PACall, ~ . -Sex)
anova(PACall2, PACall)


#~~ technician - x=38.312      df=1  p=6.031e-10

PACall3<-update(PACall, ~ . -Technician)
anova(PACall3, PACall)


#~~ MAC - x=1.7802      df=1     p=0.1821

PACall4<-update(PACall, ~ . -MAC)
anova(PACall4, PACall)


#~~ PAC - x=3.2797     df=1    p=0.07014

PACall5<-update(PACall, ~ . -PAC)
anova(PACall5, PACall)


#~~ are within vs between pac or mac effects sig?

PACallcen <- lmer(cRTL ~ logAge + Sex + Technician + MeanMAC + DevMeanMAC +
                       MeanPAC + DevMeanPAC +
                    (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                  data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACallcen)


#~~ model diagnostics

plot(PACallcen)
qqnorm(resid(PACallcen))
vif.mer(PACallcen)


#~~ save model output

PACallcenmodel<-as.data.frame((summary(PACallcen)$coefficients))
PACallcenmodel<-setDT(PACallcenmodel, keep.rownames = TRUE)[]

write.csv(PACallcenmodel, file="PACMACallcen190220.csv")

PACallcenmodel2<-as.data.frame(VarCorr(PACallcen))

write.csv(PACallcenmodel2, file="PACMACallcenranef190220.csv")


#~~ convert slopes to correlation coefficients - r=beta slope*(sd x / sd y)

#~~ between MAC

0.029710*(sd(telodata.all$MeanMAC)/sd(telodata.all$cRTL))


#~~ within PAC

-0.032146*(sd(telodata.all$DevMeanPAC)/sd(telodata.all$cRTL))


#~~ using the nakagawa and cuthill 2007 equation - t/sqrt(t2+df)

#~~ between MAC

2.737/sqrt((2.737*2.737)+2353)


#~~ within PAC

-2.831/sqrt((-2.831*-2.831)+2353)


#~~ significance of fixed effects

#~~ age - x=102.25      df=1  p<2.2e-16

PACallcen2<-update(PACallcen, ~ . -logAge)
anova(PACallcen2, PACallcen)


#~~ sex- x=0.3742      df=1     p=0.5407

PACallcen3<-update(PACallcen, ~ . -Sex)
anova(PACallcen3, PACallcen)


#~~ tech - x=35.954      df=1  p=2.021e-09

PACallcen4<-update(PACallcen, ~ . -Technician)
anova(PACallcen4, PACallcen)


#~~ mean mac - x=7.4815      df=1   p=0.006233

PACallcen5<-update(PACallcen, ~ . -MeanMAC)
anova(PACallcen5, PACallcen)


#~~ dev mean mac - x=0.7259      df=1     p=0.3942

PACallcen6<-update(PACallcen, ~ . -DevMeanMAC)
anova(PACallcen6, PACallcen)


#~~ mean pac - x=0.0105     df=1     p=0.9184

PACallcen7<-update(PACallcen, ~ . -MeanPAC)
anova(PACallcen7, PACallcen)


#~~ dev mean pac - x=8.0318      df=1   p=0.004596

PACallcen8<-update(PACallcen, ~ . -DevMeanPAC)
anova(PACallcen8, PACallcen)


#~~ graph for manuscript - figure 1 - plot within PAC and between MAC effects 

PACallcen <- lmer(cRTL ~ logAge + Sex + Technician + MeanMAC + DevMeanMAC +
                       MeanPAC + DevMeanPAC +
                       (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                     data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACallcen)


#~~ make dataframe - for within PAC effects 
#~~ all continuous variables set to mean values except for within PAC effects
#~~ estimated for males and tech '1'

predDFPAC <- data.frame (logAge = mean(telodata.all$logAge),
                         Sex = '1',
                         Technician='1',
                         MeanMAC = mean(telodata.all$MeanMAC), 
                         DevMeanMAC = mean(telodata.all$DevMeanMAC),   
                         MeanPAC = mean(telodata.all$MeanPAC),  
                         DevMeanPAC = unique(sort(telodata.all$DevMeanPAC)))


#~~ predict and add to dataframe

predPAC <- predictSE (PACallcen, predDFPAC, re.form=NA, type='response', se.fit=T)
predDFPAC$pred<-predPAC$fit
predDFPAC$upper<-predPAC$fit+predPAC$se.fit
predDFPAC$lower<-predPAC$fit-predPAC$se.fit


#~~ graph

WnPACgraph<-ggplot(telodata.all, aes(x=DevMeanPAC, y=cRTL))+ geom_point(colour="grey6", alpha=0.30) +
  xlab("Delta age father (years)") + ylab("Offspring relative telomere length") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=14),
        axis.title.y=element_text(colour="grey6", size=14),
        axis.line=element_line(colour="grey6"),
        text = element_text(size=14, colour="grey26")) +
  geom_line(data=predDFPAC, aes(x=DevMeanPAC, y=pred), colour="midnightblue", size=0.8) + 
  geom_line(data=predDFPAC, aes(x=DevMeanPAC, y=upper), colour="midnightblue", linetype="longdash", size=0.4) +
  geom_line(data=predDFPAC, aes(x=DevMeanPAC, y=lower), colour="midnightblue", linetype="longdash", size=0.4) 

WnPACgraph


#~~ make dataframe - for between MAC effects
#~~ all continuous variables set to mean values except for between MAC effects
#~~ estimated for males and tech '1'

predDFMAC <- data.frame (logAge = mean(telodata.all$logAge),
                         Sex = '1',
                         Technician='1',
                         MeanMAC = unique(sort(telodata.all$MeanMAC)), 
                         DevMeanMAC = mean(telodata.all$DevMeanMAC),   
                         MeanPAC = mean(telodata.all$MeanPAC),  
                         DevMeanPAC = mean(telodata.all$DevMeanPAC))


#~~ predict and add to dataframe

predMAC <- predictSE (PACallcen, predDFMAC, re.form=NA, type='response', se.fit=T)
predDFMAC$pred<-predMAC$fit
predDFMAC$upper<-predMAC$fit+predMAC$se.fit
predDFMAC$lower<-predMAC$fit-predMAC$se.fit

BnMACgraph<-ggplot(telodata.all, aes(x=MeanMAC, y=cRTL))+ geom_point(colour="grey6", alpha=0.30) +
  xlab("Mean maternal age (years)") + ylab("Offspring relative telomere length") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=14),
        axis.title.y=element_text(colour="grey6", size=14),
        axis.line=element_line(colour="grey6"),
        text = element_text(size=14, colour="grey26")) +
  geom_line(data=predDFMAC, aes(x=MeanMAC, y=pred), colour="midnightblue", size=0.8) + 
  geom_line(data=predDFMAC, aes(x=MeanMAC, y=upper), colour="midnightblue", linetype="longdash", size=0.4) +
  geom_line(data=predDFMAC, aes(x=MeanMAC, y=lower), colour="midnightblue", linetype="longdash", size=0.4) 

BnMACgraph

WnBnGraphs <- plot_grid(WnPACgraph, BnMACgraph, ncol=2, nrow=1, labels="AUTO")

pdf("Figure1.pdf", width=9, height=4.5)
print(WnBnGraphs)
dev.off()


#~~ are the slopes different?

PACallcenb <- lmer(cRTL ~ logAge + Sex + Technician + MeanMAC + MAC +
                       MeanPAC + PAC +
                       (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                     data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACallcenb)


#~~ model diagnostics

plot(PACallcenb)
qqnorm(resid(PACallcenb))
vif.mer(PACallcenb)


#~~ save model output

PACallcenbmodel<-as.data.frame((summary(PACallcenb)$coefficients))
PACallcenbmodel<-setDT(PACallcenbmodel, keep.rownames = TRUE)[]

write.csv(PACallcenbmodel, file="PACMACallcenb190220.csv")

PACallcenbmodel2<-as.data.frame(VarCorr(PACallcenb))

write.csv(PACallcenbmodel2, file="PACMACallcenbranef190220.csv")


#~~ slopes diff for mac - x=6.002      df=1    p=0.01429 

PACallcenb2<-update(PACallcenb, ~ . -MeanMAC)
anova(PACallcenb2, PACallcenb)


#~~ and pac - x=5.1018      df=1     p=0.0239

PACallcenb3<-update(PACallcenb, ~ . -MeanPAC)
anova(PACallcenb3, PACallcenb)


#~~ age - x=102.25      df=1  p<2.2e-16

PACallcenb4<-update(PACallcenb, ~ . -logAge)
anova(PACallcenb4, PACallcenb)


#~~ sex - x=0.3742      df=1     p=0.5407

PACallcenb5<-update(PACallcenb, ~ . -Sex)
anova(PACallcenb5, PACallcenb)


#~~ tech - x=35.954      df=1  p=2.021e-09

PACallcenb6<-update(PACallcenb, ~ . -Technician)
anova(PACallcenb6, PACallcenb)


#~~ mac - x=0.7259      df=1     p=0.3942

PACallcenb7<-update(PACallcenb, ~ . -MAC)
anova(PACallcenb7, PACallcenb)


#~~ pac - x=8.0318      df=1   p=0.004596 

PACallcenb8<-update(PACallcenb, ~ . -PAC)
anova(PACallcenb8, PACallcenb)


# # -----------------------------------------------------------------------

# PAC/MAC - chicks - first measure 

# # -----------------------------------------------------------------------


#~~ check age classes

table(telodata$AgeClass)


#~~ make chick dataset

telodata.ch<-droplevels(subset(telodata, AgeClass=="CH" & !is.na(AgeY) & !is.na(Sex) & !is.na(Technician) & !is.na(MAC) & !is.na(PAC)))


#~~ only 4 replicates - issues of model converging with ID as random effect - take first measure

telodata.ch <- arrange(telodata.ch, ID, AgeY)
telodata.ch <- telodata.ch [ ! duplicated (telodata.ch$ID), ]
anyDuplicated(telodata.ch$ID) 


#~~ z-transform sqrt RTL values

telodata.ch$cRTL <- (telodata.ch$sqrtRTL - mean(telodata.ch$sqrtRTL, na.rm=TRUE)) / sd(telodata.ch$sqrtRTL, na.rm=TRUE)


#~~ model - age has a lower aic than log age

PACchick <- lmer(cRTL ~ AgeY + Sex + Technician + MAC + PAC +
                   (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                 data=telodata.ch, control = lmerControl(optimizer ="bobyqa"))


# PACchick2 <- lmer(cRTL ~ logAge + Sex + Technician + MAC + PAC +
#                    (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
#                  data=telodata.ch, control = lmerControl(optimizer ="bobyqa"))
# 
# AIC(PACchick)
# AIC(PACchick2)

summary(PACchick)


#~~ model diagnostics 

plot(PACchick)
qqnorm(resid(PACchick)) + qqline(resid(PACchick))
vif.mer(PACchick)


#~~ save model output

PACchicksmodel<-as.data.frame((summary(PACchick)$coefficients))
PACchicksmodel<-setDT(PACchicksmodel, keep.rownames = TRUE)[]

write.csv(PACchicksmodel, file="PACMACchick190220.csv")

PACchicksmodel2<-as.data.frame(VarCorr(PACchick))

write.csv(PACchicksmodel2, file="PACMACchicksranef190220.csv")


#~~ check sig of fixed effects

#~~ age - x=5.8118      df=1    p=0.01592

PACchick1<-update(PACchick, ~ . -AgeY)
anova(PACchick1, PACchick)


#~~ sex - x=0.4004      df=1     p=0.5269

PACchick2<-update(PACchick, ~ . -Sex)
anova(PACchick2, PACchick)


#~~ technician - x=2.961      df=1    p=0.08529

PACchick3<-update(PACchick, ~ . -Technician)
anova(PACchick3, PACchick)


#~~ MAC - x=1e-04      df=1     p=0.9922

PACchick4<-update(PACchick, ~ . -MAC)
anova(PACchick4, PACchick)


#~~ PAC - x=0.1641      df=1     p=0.6854

PACchick5<-update(PACchick, ~ . -PAC)
anova(PACchick5, PACchick)



# # -----------------------------------------------------------------------

# PAC/MAC - juveniles

# # -----------------------------------------------------------------------


#~~ make juvenile dataset (<1 year)

telodata.juv<-droplevels(subset(telodata, AgeY<1 & !is.na(AgeY) & !is.na(Sex) & !is.na(Technician) & !is.na(MAC) & !is.na(PAC)))


#~~ z-transform sqrt RTL values

telodata.juv$cRTL <- (telodata.juv$sqrtRTL - mean(telodata.juv$sqrtRTL, na.rm=TRUE)) / sd(telodata.juv$sqrtRTL, na.rm=TRUE)


#~~ logAge has lower AIC than age - model would not converge with dad as a random effect

PACjuv <- lmer(cRTL ~ logAge + Sex + Technician + MAC + PAC +
                   (1|ID) + (1|mum) + (1|FPID) + (1|U_PlateID),
                 data=telodata.juv, control = lmerControl(optimizer ="bobyqa"))


# PACjuv2 <- lmer(cRTL ~ AgeY + Sex + Technician + MAC + PAC +
#                  (1|ID) + (1|mum) + (1|FPID) + (1|U_PlateID),
#                data=telodata.juv, control = lmerControl(optimizer ="bobyqa"))
# 
# AIC(PACjuv)
# AIC(PACjuv2)

summary(PACjuv)


#~~ model diagnostics

plot(PACjuv)
qqnorm(resid(PACjuv)) + qqline(resid(PACjuv))
vif.mer(PACjuv)


#~~ save model outputs 

PACjuvmodel<-as.data.frame((summary(PACjuv)$coefficients))
PACjuvmodel<-setDT(PACjuvmodel, keep.rownames = TRUE)[]

write.csv(PACjuvmodel, file="PACMACjuv190220.csv")

PACjuvmodel2<-as.data.frame(VarCorr(PACjuv))

write.csv(PACjuvmodel2, file="PACMACjuvranef190220.csv")


#~~ check sig of fixed effects

#~~ age - x=66.225      df=1  p=4.023e-16

PACjuv1<-update(PACjuv, ~ . -logAge)
anova(PACjuv1, PACjuv)


#~~ sex - x=1.9327      df=1     p=0.1645

PACjuv2<-update(PACjuv, ~ . -Sex)
anova(PACjuv2, PACjuv)


#~~ technician - x=12.809      df=1  p=0.0003449

PACjuv3<-update(PACjuv, ~ . -Technician)
anova(PACjuv3, PACjuv)


#~~ MAC - x=0.1667      df=1     p=0.6831

PACjuv4<-update(PACjuv, ~ . -MAC)
anova(PACjuv4, PACjuv)


#~~ PAC - x= 0.7151      df=1     p=0.3978

PACjuv5<-update(PACjuv, ~ . -PAC)
anova(PACjuv5, PACjuv)


#~~ graph for supplementary material - figure s8 - offspring RTL (all ages/chicks/juveniles) vs PAC/MAC

MACall<-ggplot(telodata.all,aes(x=MAC,y=cRTL)) + geom_point(shape=19, alpha=0.2, colour="grey6")+
  labs(x="", y="") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="white", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) + ggtitle(expression("All ages")) +
  theme(plot.title=element_text(size=13, color="grey6")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))

PACall<-ggplot(telodata.all,aes(x=PAC,y=cRTL)) + geom_point(shape=19, alpha=0.2, colour="grey6")+
  labs(x="", y="") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="white", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) + ggtitle(expression("/")) +
  theme(plot.title=element_text(size=13, color="white")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))

MACch<-ggplot(telodata.ch,aes(x=MAC,y=cRTL)) + geom_point(shape=19, alpha=0.4, colour="grey6")+
  labs(x="", y="Offspring relative telomere length") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="grey6", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) + ggtitle(expression("Chicks")) +
  theme(plot.title=element_text(size=13, color="grey6")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))

PACch<-ggplot(telodata.ch,aes(x=PAC,y=cRTL)) + geom_point(shape=19, alpha=0.4, colour="grey6")+
  labs(x="", y="") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="white", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) + ggtitle(expression("/")) +
  theme(plot.title=element_text(size=13, color="white")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))

MACjuv<-ggplot(telodata.juv,aes(x=MAC,y=cRTL)) + geom_point(shape=19, alpha=0.4, colour="grey6")+
  labs(x="Maternal age at conception", y="") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="white", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(colour="grey16")) + ggtitle(expression("Juveniles")) +
  theme(plot.title=element_text(size=13, color="grey6")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))

PACjuv<-ggplot(telodata.juv,aes(x=PAC,y=cRTL)) + geom_point(shape=19, alpha=0.4, colour="grey6")+
  labs(x="Paternal age at conception", y="") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) + 
  theme(axis.title.y=element_text(colour="white", size=13, vjust=3)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) + ggtitle(expression("/")) +
  theme(plot.title=element_text(size=13, color="white")) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15))


plot_grid(MACall, PACall, MACch, PACch, MACjuv, PACjuv, ncol=2, nrow=3)


MACPACgraphs<-plot_grid(MACall, PACall, MACch, PACch, MACjuv, PACjuv, ncol=2, nrow=3, align="hv", labels="AUTO", hjust=-1.8)

tiff("FigureS8_PACMACGraphs.tif", width=3000, height=3500, res=400)

MACPACgraphs

dev.off()


# # -----------------------------------------------------------------------

# RESULTS WITH UNTRANSFORMED DATA - ALL AGES 

# # -----------------------------------------------------------------------


#~~ graph for supplementary material - figures s4 - histograms of transformed vs untransformed RTL

rtldistgraph<-ggplot(telodata.all, aes(x=RTL)) +
  geom_histogram(alpha=.5, position="identity", fill="grey40") +
  xlab("RTL") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey16", size=16)) + 
  theme(axis.title.y=element_text(colour="grey16", size=16)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=17, colour="grey26")) 

crtldistgraph<-ggplot(telodata.all, aes(x=cRTL)) +
  geom_histogram(alpha=.5, position="identity", fill="grey40") +
  xlab("Square root and z-transformed RTL") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey16", size=16)) + 
  theme(axis.title.y=element_text(colour="grey16", size=16)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=17, colour="grey26")) 


plot_grid(rtldistgraph, crtldistgraph, ncol = 2, nrow=1)

hist.figs<-plot_grid(rtldistgraph, crtldistgraph, ncol = 2, nrow=1, align="hv", labels="AUTO")

tiff("FigureS4_RTLcRTLhistograms.tif", width=4500, height=2200, res=400)

hist.figs

dev.off()


#~~ model

PACall.nt <- lmer(RTL ~ logAge + Sex + Technician + MAC + PAC +
                    (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                  data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACall.nt)


#~~ model diagnostics 

plot(PACall.nt)
qqnorm(resid(PACall.nt))
vif.mer(PACall.nt)


#~~ save model output

PACallmodel.nt<-as.data.frame((summary(PACall.nt)$coefficients))
PACallmodel.nt<-setDT(PACallmodel.nt, keep.rownames = TRUE)[]

write.csv(PACallmodel.nt, file="PACMACall_nont_111120.csv")

PACallmodel2.nt<-as.data.frame(VarCorr(PACall.nt))

write.csv(PACallmodel2.nt, file="PACMACallranef_nont_111120.csv")


#~~ check sig of fixed effects

#~~ age - x=93.926     df=1  p<2.2e-16

PACall.nt1<-update(PACall.nt, ~ . -logAge)
anova(PACall.nt1, PACall.nt)


#~~ sex - x=0.1241      df=1     p=0.7246

PACall.nt2<-update(PACall.nt, ~ . -Sex)
anova(PACall.nt2, PACall.nt)


#~~ technician - x=32.263      df=1  p=1.346e-08

PACall.nt3<-update(PACall.nt, ~ . -Technician)
anova(PACall.nt3, PACall.nt)


#~~ MAC - x=1.7393     df=1     p=0.1872

PACall.nt4<-update(PACall.nt, ~ . -MAC)
anova(PACall.nt4, PACall.nt)


#~~ PAC - x=3.5898      df=1    p=0.05814

PACall.nt5<-update(PACall.nt, ~ . -PAC)
anova(PACall.nt5, PACall.nt)


#~~ are within vs between pac or mac effects sig?

PACallcen.nt <- lmer(RTL ~ logAge + Sex + Technician + MeanMAC + DevMeanMAC +
                       MeanPAC + DevMeanPAC +
                       (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                     data=telodata.all, control = lmerControl(optimizer ="bobyqa"))

summary(PACallcen.nt)


#~~ model diagnostics 

plot(PACallcen.nt)
qqnorm(resid(PACallcen.nt))
vif.mer(PACallcen.nt)


#~~ save model output

PACallcenmodel.nt<-as.data.frame((summary(PACallcen.nt)$coefficients))
PACallcensmodel.nt<-setDT(PACallcenmodel.nt, keep.rownames = TRUE)[]

write.csv(PACallcensmodel.nt, file="PACMACallcens_nont_111120.csv")

PACallcenmodel.nt2<-as.data.frame(VarCorr(PACallcen.nt))

write.csv(PACallcenmodel.nt2, file="PACMACallcensranef_nont_111120.csv")


#~~ sig of fixed effects

#~~ age - x=101.02      df=1  p<2.2e-16

PACallcen.nt2<-update(PACallcen.nt, ~ . -logAge)
anova(PACallcen.nt2, PACallcen.nt)


#~~ sex - x=0.1637      df=1     p=0.6857

PACallcen.nt3<-update(PACallcen.nt, ~ . -Sex)
anova(PACallcen.nt3, PACallcen.nt)


#~~ tech - x=30.196      df=1  p=3.905e-08

PACallcen.nt4<-update(PACallcen.nt, ~ . -Technician)
anova(PACallcen.nt4, PACallcen.nt)


#~~ mean mac - x=6.8228      df=1     p=0.009

PACallcen.nt5<-update(PACallcen.nt, ~ . -MeanMAC)
anova(PACallcen.nt5, PACallcen.nt)


#~~ dev mean mac - x=0.5843      df=1     p=0.4446

PACallcen.nt6<-update(PACallcen.nt, ~ . -DevMeanMAC)
anova(PACallcen.nt6, PACallcen.nt)


#~~ mean pac - x=0      df=1    p=0.9961

PACallcen.nt7<-update(PACallcen.nt, ~ . -MeanPAC)
anova(PACallcen.nt7, PACallcen.nt)


#~~ dev mean pac - x=8.4085      df=1  p=0.003735

PACallcen.nt8<-update(PACallcen.nt, ~ . -DevMeanPAC)
anova(PACallcen.nt8, PACallcen.nt)


#~~ are the slopes different? 

PACallcenb.nt <- lmer(RTL ~ logAge + Sex + Technician + MeanMAC + MAC +
                        MeanPAC + PAC +
                        (1|ID) + (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                      data=telodata.all, control = lmerControl(optimizer ="bobyqa"))


summary(PACallcenb.nt)


#~~ model diagnostics

plot(PACallcenb.nt)
qqnorm(resid(PACallcenb.nt))
vif.mer(PACallcenb.nt)


#~~ save model output

PACallcenb.nt.model<-as.data.frame((summary(PACallcenb.nt)$coefficients))
PACallcenb.nt.model<-setDT(PACallcenb.nt.model, keep.rownames = TRUE)[]

write.csv(PACallcenb.nt.model, file="PACMACallcenb_nont_141120.csv")

PACallcenb.nt.model2<-as.data.frame(VarCorr(PACallcenb.nt))

write.csv(PACallcenb.nt.model2, file="PACMACallcenbranef_nont_141120.csv")


#~~ slopes diff for mac - x=5.2817      df=1    p=0.02155

PACallcenb.nt2<-update(PACallcenb.nt, ~ . -MeanMAC)
anova(PACallcenb.nt2, PACallcenb.nt)


#~~ and pac - x=5.1422      df=1    p=0.02335

PACallcenb.nt3<-update(PACallcenb.nt, ~ . -MeanPAC)
anova(PACallcenb.nt3, PACallcenb.nt)


#~~ age - x=101.02      df=1  p<2.2e-16

PACallcenb.nt4<-update(PACallcenb.nt, ~ . -logAge)
anova(PACallcenb.nt4, PACallcenb.nt)


#~~ sex - x=0.1637     df=1     p=0.6857

PACallcenb.nt5<-update(PACallcenb.nt, ~ . -Sex)
anova(PACallcenb.nt5, PACallcenb.nt)


#~~ tech -x=30.196      df=1  p=3.905e-08

PACallcenb.nt6<-update(PACallcenb.nt, ~ . -Technician)
anova(PACallcenb.nt6, PACallcenb.nt)


#~~ mac - x=0.5843      df=1     p=0.4446

PACallcenb.nt7<-update(PACallcenb.nt, ~ . -MAC)
anova(PACallcenb.nt7, PACallcenb.nt)


#~~ pac - x=8.4085      df=1   p=0.003735

PACallcenb.nt8<-update(PACallcenb.nt, ~ . -PAC)
anova(PACallcenb.nt8, PACallcenb.nt)


# # -----------------------------------------------------------------------

# RESULTS WITH UNTRANSFORMED DATA - CHICKS

# # -----------------------------------------------------------------------


#~~ model 

PACchick.nt <- lmer(RTL ~ AgeY + Sex + Technician + MAC + PAC +
                   (1|mum) + (1|dad) + (1|FPID) + (1|U_PlateID),
                 data=telodata.ch, control = lmerControl(optimizer ="bobyqa"))


summary(PACchick.nt)


#~~ model diagnostics 

plot(PACchick.nt)
qqnorm(resid(PACchick.nt)) + qqline(resid(PACchick.nt))
vif.mer(PACchick.nt)


#~~ save model output

PACchicksmodel.nt<-as.data.frame((summary(PACchick.nt)$coefficients))
PACchicksmodel.nt<-setDT(PACchicksmodel.nt, keep.rownames = TRUE)[]

write.csv(PACchicksmodel.nt, file="PACMACchick_nont_141120.csv")

PACchicksmodel.nt2<-as.data.frame(VarCorr(PACchick.nt))

write.csv(PACchicksmodel.nt2, file="PACMACchicksranef_nont_141120.csv")


#~~ check sig of fixed effects

#~~ age - x=6.787     df=1   p=0.009183

PACchick.nt1<-update(PACchick.nt, ~ . -AgeY)
anova(PACchick.nt1, PACchick.nt)


#~~ sex - x=0.3928      df=1     p=0.5308

PACchick.nt2<-update(PACchick.nt, ~ . -Sex)
anova(PACchick.nt2, PACchick.nt)


#~~ technician - x=2.5977      df=1      p=0.107

PACchick.nt3<-update(PACchick.nt, ~ . -Technician)
anova(PACchick.nt3, PACchick.nt)


#~~ MAC - x=0.0025      df=1       p=0.96

PACchick.nt4<-update(PACchick.nt, ~ . -MAC)
anova(PACchick.nt4, PACchick.nt)


#~~ PAC - x=0.0246      df=1    p=0.8754

PACchick.nt5<-update(PACchick.nt, ~ . -PAC)
anova(PACchick.nt5, PACchick.nt)


# # -----------------------------------------------------------------------

# RESULTS WITH UNTRANSFORMED DATA - JUVENILES

# # -----------------------------------------------------------------------


#~~ model 

PACjuv.nt <- lmer(RTL ~ logAge + Sex + Technician + MAC + PAC +
                 (1|ID) + (1|mum) + (1|FPID) + (1|U_PlateID),
               data=telodata.juv, control = lmerControl(optimizer ="bobyqa"))

summary(PACjuv.nt)


#~~ model diagnostics

plot(PACjuv.nt)
qqnorm(resid(PACjuv.nt)) + qqline(resid(PACjuv.nt))
vif.mer(PACjuv.nt)


#~~ save model output 

PACjuvmodel.nt<-as.data.frame((summary(PACjuv.nt)$coefficients))
PACjuvmodel.nt<-setDT(PACjuvmodel.nt, keep.rownames = TRUE)[]

write.csv(PACjuvmodel.nt, file="PACMACjuv_nont_141120.csv")


PACjuvmodel.nt2<-as.data.frame(VarCorr(PACjuv.nt))

write.csv(PACjuvmodel.nt2, file="PACMACjuvranef_nont_141120.csv")


#~~ check sig of fixed effects

#~~ age - x=64.221      df=1  p=1.112e-15

PACjuv.nt1<-update(PACjuv.nt, ~ . -logAge)
anova(PACjuv.nt1, PACjuv.nt)


#~~ sex - x=1.0408      df=1     p=0.3076

PACjuv.nt2<-update(PACjuv.nt, ~ . -Sex)
anova(PACjuv.nt2, PACjuv.nt)


#~~ technician - x=12.071      df=1  p=0.0005121

PACjuv.nt3<-update(PACjuv.nt, ~ . -Technician)
anova(PACjuv.nt3, PACjuv.nt)


#~~ MAC - x=0.1624      df=1      p=0.687

PACjuv.nt4<-update(PACjuv.nt, ~ . -MAC)
anova(PACjuv.nt4, PACjuv.nt)


#~~ PAC - x=0.914      df=1     p=0.339

PACjuv.nt5<-update(PACjuv.nt, ~ . -PAC)
anova(PACjuv.nt5, PACjuv.nt)

