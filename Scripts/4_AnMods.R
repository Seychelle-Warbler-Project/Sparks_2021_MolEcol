# # -----------------------------------------------------------------------

# ----- Telomere heritability manuscript
# ----- AMS
# ----- Animal models
# ----- The models

# # -----------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(MCMCglmm)
library(MasterBayes)
library(dplyr)
library(pedantics)
library(asreml)
library(nadiv)
library(ggplot2)
library(merTools)


#~~ data 

telodata <- read.csv("Sparks_Telo_h2_MS_Data_2020.csv", header=T)

colnames(telodata)[which(names(telodata) == "BirdID")] <- "id"

telodata$id <- as.factor(telodata$id)
telodata$animal <- telodata$id
telodata$Sex <- as.factor(telodata$Sex)
telodata$Technician <- as.factor(telodata$Technician)
telodata$FPID <- as.factor(telodata$FPID)
telodata$Terr <- as.factor(telodata$Terr)
telodata$mum <- as.factor(telodata$mum)
telodata$dad <- as.factor(telodata$dad)
telodata$BrF <- as.factor(telodata$BrF)
telodata$BrM <- as.factor(telodata$BrM)
telodata$U_PlateID <- as.factor(telodata$U_PlateID)
telodata$BirthFPID <- as.factor(telodata$BirthFPID)
telodata$sqrtRTL <- sqrt(telodata$RTL)
telodata$logAge <- log10(telodata$AgeY)


#~~ ped

ped <- read.csv("Telo_h2_Pruned_Pedigree_2020.csv", header=T)


#~~  list dams/sires before offspring

ped <- orderPed(ped) 


#~~ subset telomere dataset to only individuals in ped for anmods to run

ped$id<-as.factor(ped$id)
ped$mum<-as.factor(ped$mum)
ped$dad<-as.factor(ped$dad)

telodata.ped<-telodata[telodata$id%in%ped$id,]
telodata.ped<-droplevels(telodata.ped)


#~~ to clean up workspace 

KeepThese <- c("Results", "Results2", "Results3", "Results4", "Results5", "ped", "telodata.ped",
               "Results6", "Results7", "Results8", "Results9", "Results10", "KeepThese", "prior2")


# # -----------------------------------------------------------------------

# Model 1

# # -----------------------------------------------------------------------


prior1 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
                R=list(V=1, n=0.002))

model1<-MCMCglmm(RTL~1, random=~id, data=telodata.ped,
                 nitt=1200000, thin=500, burnin=200000, prior=prior1)

save(model1, file="AnMod_MCMCglmm_1_280220.RData")

load(file="AnMod_MCMCglmm_1_280220.RData")


#~~ diagnostics

plot(model1$Sol) 
plot(model1$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model1$Sol)
autocorr(model1$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model1$Sol)
effectiveSize(model1$VCV)


#~~ heidel diagnostics

heidel.diag(model1$Sol)
heidel.diag(model1$VCV)


#~~ geweke diagnostics

geweke.diag(model1$VCV)
geweke.diag(model1$Sol)


#~~ for model output

posterior.mode(model1$Sol)
HPDinterval(model1$Sol, prob=0.95)
posterior.mode(model1$VCV)
HPDinterval(model1$VCV, prob=0.95)


#~~ variance estimates

idvar1 <- posterior.mode(model1$VCV[,"id"])
idvarcri1 <- HPDinterval(model1$VCV[,"id"], 0.95)

resvar1 <- posterior.mode(model1$VCV[,"units"])
resvarcri1 <- HPDinterval(model1$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

rep1 <- model1$VCV[,"id"]/(model1$VCV[,"id"]+model1$VCV[,"units"])
repmodel1 <- posterior.mode(rep1)
repmodel1cri <- HPDinterval(rep1, 0.95)

repmodel1
repmodel1cri

res1 <- model1$VCV[,"units"]/(model1$VCV[,"id"]+model1$VCV[,"units"])
resmodel1 <- posterior.mode(res1)
resmodel1cri <- HPDinterval(res1, 0.95)


#~~ for results

Variance <- c("ID", "Residual")
VarEst <- rbind(idvar1, resvar1)
VarCri <- rbind(idvarcri1, resvarcri1)
PropEst <- rbind(repmodel1, resmodel1)
PropCri <- rbind(repmodel1cri, resmodel1cri)

Results <- as.data.frame(cbind(Variance, VarEst, VarCri, PropEst, PropCri))
Results$Model <- "1"


# # -----------------------------------------------------------------------

# Model 2

# # -----------------------------------------------------------------------


prior2 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                       G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
                R=list(V=1, n=0.002))

model2<-MCMCglmm(RTL~1, random=~animal + id, pedigree=ped, data=telodata.ped,
                 nitt=1200000, thin=500, burnin=200000, prior=prior2)

save(model2, file="AnMod_MCMCglmm_2_260220v2.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_2_260220v2.RData")


#~~ diagnostics

plot(model2$Sol) 
plot(model2$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model2$Sol)
autocorr(model2$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model2$Sol)
effectiveSize(model2$VCV)


#~~ heidel

heidel.diag(model2$Sol)
heidel.diag(model2$VCV)


#~~ geweke

geweke.diag(model2$VCV)
geweke.diag(model2$Sol)


#~~ for model output

posterior.mode(model2$Sol)
HPDinterval(model2$Sol, prob=0.95)
posterior.mode(model2$VCV)
HPDinterval(model2$VCV, prob=0.95)


#~~ variance estimates

addgenvar2 <- posterior.mode(model2$VCV[,"animal"])
addgencri2 <- HPDinterval(model2$VCV[,"animal"], 0.95)

perenvvar2 <- posterior.mode(model2$VCV[,"id"])
perenvcri2 <- HPDinterval(model2$VCV[,"id"], 0.95)

resvar2 <- posterior.mode(model2$VCV[,"units"])
resvarcri2 <- HPDinterval(model2$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen2 <- model2$VCV[,"animal"]/(model2$VCV[,"id"]+model2$VCV[,"animal"]+model2$VCV[,"units"])
addgenmodel2 <- posterior.mode(addgen2)
addgenmodel2cri <- HPDinterval(addgen2, 0.95)

addgenmodel2
addgenmodel2cri

perenv2 <- model2$VCV[,"id"]/(model2$VCV[,"id"]+model2$VCV[,"animal"]+model2$VCV[,"units"])
perenvmodel2 <- posterior.mode(perenv2)
perenvmodel2cri <- HPDinterval(perenv2, 0.95)

res2 <- model2$VCV[,"units"]/(model2$VCV[,"id"]+model2$VCV[,"animal"]+model2$VCV[,"units"])
resmodel2 <- posterior.mode(res2)
resmodel2cri <- HPDinterval(res2, 0.95)


#~~  repeatability 

rep2 <- (model2$VCV[,"id"]+model2$VCV[,"animal"])/(model2$VCV[,"id"]+model2$VCV[,"animal"]+model2$VCV[,"units"])
repmodel2 <- posterior.mode(rep2)
repmodel2 <- HPDinterval(rep2, 0.95)


#~~  for results

Variance <- c("AddGen", "PerEnv", "Residual")
VarEst2 <- rbind(addgenvar2, perenvvar2, resvar2)
VarCri2 <- rbind(addgencri2, perenvcri2, resvarcri2)
PropEst2 <- rbind(addgenmodel2, perenvmodel2, resmodel2)
PropCri2 <- rbind(addgenmodel2cri, perenvmodel2cri, resmodel2cri)

Results2 <- as.data.frame(cbind(Variance, VarEst2, VarCri2, PropEst2, PropCri2))
Results2$Model <- "2"


# # -----------------------------------------------------------------------

# Model 3

# # -----------------------------------------------------------------------


model3<-MCMCglmm(RTL~Sex + logAge + Technician, random=~animal + id, pedigree=ped, data=telodata.ped,
                 nitt=1200000, thin=500, burnin=200000, prior=prior2)

save(model3, file="AnMod_MCMCglmm_3_270220.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_3_270220.RData")


#~~ diagnostics

plot(model3$Sol) 
plot(model3$VCV)


#~~ correlation among thinned samples must be low (<0.05) 

autocorr(model3$Sol)
autocorr(model3$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model3$Sol)
effectiveSize(model3$VCV)


#~~ heidel

heidel.diag(model3$Sol)
heidel.diag(model3$VCV)


#~~ geweke

geweke.diag(model3$VCV)
geweke.diag(model3$Sol)


#~~ variance estimates

addgenvar3 <- posterior.mode(model3$VCV[,"animal"])
addgencri3 <- HPDinterval(model3$VCV[,"animal"], 0.95)

perenvvar3 <- posterior.mode(model3$VCV[,"id"])
perenvcri3 <- HPDinterval(model3$VCV[,"id"], 0.95)

resvar3 <- posterior.mode(model3$VCV[,"units"])
resvarcri3 <- HPDinterval(model3$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen3 <- model3$VCV[,"animal"]/(model3$VCV[,"id"]+model3$VCV[,"animal"]+model3$VCV[,"units"])
addgenmodel3 <- posterior.mode(addgen3)
addgenmodel3cri <- HPDinterval(addgen3, 0.95)

perenv3 <- model3$VCV[,"id"]/(model3$VCV[,"id"]+model3$VCV[,"animal"]+model3$VCV[,"units"])
perenvmodel3 <- posterior.mode(perenv3)
perenvmodel3cri <- HPDinterval(perenv3, 0.95)

res3 <- model3$VCV[,"units"]/(model3$VCV[,"id"]+model3$VCV[,"animal"]+model3$VCV[,"units"])
resmodel3 <- posterior.mode(res3)
resmodel3cri <- HPDinterval(res3, 0.95)


#~~ repeatability

rep3 <- (model3$VCV[,"id"]+model3$VCV[,"animal"])/(model3$VCV[,"id"]+model3$VCV[,"animal"]+model3$VCV[,"units"])
repmodel3 <- posterior.mode(rep3)
repmodel3cri <- HPDinterval(rep3, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Residual")
VarEst3 <- rbind(addgenvar3, perenvvar3, resvar3)
VarCri3 <- rbind(addgencri3, perenvcri3, resvarcri3)
PropEst3 <- rbind(addgenmodel3, perenvmodel3, resmodel3)
PropCri3 <- rbind(addgenmodel3cri, perenvmodel3cri, resmodel3cri)

Results3 <- as.data.frame(cbind(Variance, VarEst3, VarCri3, PropEst3, PropCri3))
Results3$Model <- "3"


# # -----------------------------------------------------------------------

# Model 4

# # -----------------------------------------------------------------------


prior3 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), 
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model4<-MCMCglmm(RTL~Sex + logAge + Technician, random=~animal + id + U_PlateID, pedigree=ped, data=telodata.ped,
                 nitt=2400000, thin=1000, burnin=400000, prior=prior3)

save(model4, file="AnMod_MCMCglmm_4a_230320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_4a_230320.RData")


#~~ diagnostics

plot(model4$Sol) 
plot(model4$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model4$Sol)
autocorr(model4$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model4$Sol)
effectiveSize(model4$VCV)


#~~ heidel

heidel.diag(model4$Sol)
heidel.diag(model4$VCV)


#~~ geweke

geweke.diag(model4$VCV)
geweke.diag(model4$Sol)


#~~ variance estimates

addgenvar4 <- posterior.mode(model4$VCV[,"animal"])
addgencri4 <- HPDinterval(model4$VCV[,"animal"], 0.95)

perenvvar4 <- posterior.mode(model4$VCV[,"id"])
perenvcri4 <- HPDinterval(model4$VCV[,"id"], 0.95)

platevar4 <- posterior.mode(model4$VCV[,"U_PlateID"])
platecri4 <- HPDinterval(model4$VCV[,"U_PlateID"], 0.95)

resvar4 <- posterior.mode(model4$VCV[,"units"])
resvarcri4 <- HPDinterval(model4$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen4 <- model4$VCV[,"animal"]/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                    model4$VCV[,"U_PlateID"]+model4$VCV[,"units"])
addgenmodel4 <- posterior.mode(addgen4)
addgenmodel4cri <- HPDinterval(addgen4, 0.95)

perenv4 <- model4$VCV[,"id"]/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                model4$VCV[,"U_PlateID"]+model4$VCV[,"units"])
perenvmodel4 <- posterior.mode(perenv4)
perenvmodel4cri <- HPDinterval(perenv4, 0.95)

plate4 <- model4$VCV[,"U_PlateID"]/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                model4$VCV[,"U_PlateID"]+model4$VCV[,"units"])
platemodel4 <- posterior.mode(plate4)
platemodel4cri <- HPDinterval(plate4, 0.95)

res4 <- model4$VCV[,"units"]/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                model4$VCV[,"U_PlateID"]+model4$VCV[,"units"])
resmodel4 <- posterior.mode(res4)
resmodel4cri <- HPDinterval(res4, 0.95)


#~~ repeatability

rep4 <- (model4$VCV[,"id"]+model4$VCV[,"animal"])/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                                     model4$VCV[,"U_PlateID"]+model4$VCV[,"units"])
repmodel4 <- posterior.mode(rep4)
repmodel4cri <- HPDinterval(rep4, 0.95)


#~~ h2 without plate

addgen4woplate <- model4$VCV[,"animal"]/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                    model4$VCV[,"units"])
addgenmodel4woplate <- posterior.mode(addgen4woplate)
addgenmodel4criwoplate <- HPDinterval(addgen4woplate, 0.95)


#~~ rep without plate

rep4woplate <- (model4$VCV[,"id"]+model4$VCV[,"animal"])/(model4$VCV[,"id"]+model4$VCV[,"animal"]+
                                                     model4$VCV[,"units"])
repmodel4woplate <- posterior.mode(rep4woplate)
repmodel4criwoplate <- HPDinterval(rep4woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Residual")
VarEst4 <- rbind(addgenvar4, perenvvar4, platevar4, resvar4)
VarCri4 <- rbind(addgencri4, perenvcri4, platecri4, resvarcri4)
PropEst4 <- rbind(addgenmodel4, perenvmodel4, platemodel4, resmodel4)
PropCri4 <- rbind(addgenmodel4cri, perenvmodel4cri, platemodel4cri, resmodel4cri)

Results4 <- as.data.frame(cbind(Variance, VarEst4, VarCri4, PropEst4, PropCri4))
Results4$Model <- "4"


# # -----------------------------------------------------------------------

# Model 5

# # -----------------------------------------------------------------------


prior4 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), 
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model5<-MCMCglmm(RTL~Sex + logAge + Technician, random=~animal + id + U_PlateID + mum, pedigree=ped, data=telodata.ped,
                 nitt=2400000, thin=1000, burnin=400000, prior=prior4)

save(model5, file="AnMod_MCMCglmm_5_040320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_5_040320.RData")


#~~ diagnostics

plot(model5$Sol) 
plot(model5$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model5$Sol)
autocorr(model5$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model5$Sol)
effectiveSize(model5$VCV)


#~~ heidel

heidel.diag(model5$Sol)
heidel.diag(model5$VCV)


#~~ geweke

geweke.diag(model5$VCV)
geweke.diag(model5$Sol)


#~~ variance estimates

addgenvar5 <- posterior.mode(model5$VCV[,"animal"])
addgencri5 <- HPDinterval(model5$VCV[,"animal"], 0.95)

perenvvar5 <- posterior.mode(model5$VCV[,"id"])
perenvcri5 <- HPDinterval(model5$VCV[,"id"], 0.95)

platevar5 <- posterior.mode(model5$VCV[,"U_PlateID"])
platecri5 <- HPDinterval(model5$VCV[,"U_PlateID"], 0.95)

mumvar5 <- posterior.mode(model5$VCV[,"mum"])
mumcri5 <- HPDinterval(model5$VCV[,"mum"], 0.95)

resvar5 <- posterior.mode(model5$VCV[,"units"])
resvarcri5 <- HPDinterval(model5$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen5 <- model5$VCV[,"animal"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
addgenmodel5 <- posterior.mode(addgen5)
addgenmodel5cri <- HPDinterval(addgen5, 0.95)

perenv5 <- model5$VCV[,"id"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
perenvmodel5 <- posterior.mode(perenv5)
perenvmodel5cri <- HPDinterval(perenv5, 0.95)

plate5 <- model5$VCV[,"U_PlateID"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
platemodel5 <- posterior.mode(plate5)
platemodel5cri <- HPDinterval(plate5, 0.95)

mum5 <- model5$VCV[,"mum"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
mummodel5 <- posterior.mode(mum5)
mummodel5cri <- HPDinterval(mum5, 0.95)

res5 <- model5$VCV[,"units"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
resmodel5 <- posterior.mode(res5)
resmodel5cri <- HPDinterval(res5, 0.95)


#~~ repeatability

rep5 <- (model5$VCV[,"animal"]+model5$VCV[,"id"])/
           (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"U_PlateID"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
repmodel5 <- posterior.mode(rep5)
repmodel5cri <- HPDinterval(rep5, 0.95)


#~~ h2 wo plate

addgen5woplate <- model5$VCV[,"animal"]/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
addgenmodel5woplate <- posterior.mode(addgen5woplate)
addgenmodel5criwoplate <- HPDinterval(addgen5woplate, 0.95)


#~~ repeatability wo plate

rep5woplate <- (model5$VCV[,"animal"]+model5$VCV[,"id"])/
  (model5$VCV[,"animal"]+model5$VCV[,"id"]+model5$VCV[,"mum"]+model5$VCV[,"units"])
repmodel5woplate <- posterior.mode(rep5woplate)
repmodel5criwoplate <- HPDinterval(rep5woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Residual")
VarEst5 <- rbind(addgenvar5, perenvvar5, platevar5, mumvar5, resvar5)
VarCri5 <- rbind(addgencri5, perenvcri5, platecri5, mumcri5, resvarcri5)
PropEst5 <- rbind(addgenmodel5, perenvmodel5, platemodel5, mummodel5, resmodel5)
PropCri5 <- rbind(addgenmodel5cri, perenvmodel5cri, platemodel5cri, mummodel5cri, resmodel5cri)

Results5 <- as.data.frame(cbind(Variance, VarEst5, VarCri5, PropEst5, PropCri5))
Results5$Model <- "5"


# # -----------------------------------------------------------------------

# Model 6

# # -----------------------------------------------------------------------


prior5 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G5=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))


model6<-MCMCglmm(RTL~Sex + logAge + Technician,
                 random=~animal + id + U_PlateID + mum + dad, pedigree=ped, data=telodata.ped,
                 nitt=3600000, thin=1500, burnin=600000, prior=prior5)

save(model6, file="AnMod_MCMCglmm_6_060320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_6_060320.RData")


#~~ diagnostics

plot(model6$Sol) 
plot(model6$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model6$Sol)
autocorr(model6$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model6$Sol)
effectiveSize(model6$VCV)


#~~ heidel

heidel.diag(model6$Sol)
heidel.diag(model6$VCV)


#~~ geweke 

geweke.diag(model6$VCV)
geweke.diag(model6$Sol)


#~~ variance estimates 

addgenvar6 <- posterior.mode(model6$VCV[,"animal"])
addgencri6 <- HPDinterval(model6$VCV[,"animal"], 0.95)

perenvvar6 <- posterior.mode(model6$VCV[,"id"])
perenvcri6 <- HPDinterval(model6$VCV[,"id"], 0.95)

platevar6 <- posterior.mode(model6$VCV[,"U_PlateID"])
platecri6 <- HPDinterval(model6$VCV[,"U_PlateID"], 0.95)

mumvar6 <- posterior.mode(model6$VCV[,"mum"])
mumcri6 <- HPDinterval(model6$VCV[,"mum"], 0.95)

dadvar6 <- posterior.mode(model6$VCV[,"dad"])
dadcri6 <- HPDinterval(model6$VCV[,"dad"], 0.95)

resvar6 <- posterior.mode(model6$VCV[,"units"])
resvarcri6 <- HPDinterval(model6$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen6 <- model6$VCV[,"animal"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
addgenmodel6 <- posterior.mode(addgen6)
addgenmodel6cri <- HPDinterval(addgen6, 0.95)

perenv6 <- model6$VCV[,"id"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
perenvmodel6 <- posterior.mode(perenv6)
perenvmodel6cri <- HPDinterval(perenv6, 0.95)

plate6 <- model6$VCV[,"U_PlateID"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
platemodel6 <- posterior.mode(plate6)
platemodel6cri <- HPDinterval(plate6, 0.95)

mum6 <- model6$VCV[,"mum"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
mummodel6 <- posterior.mode(mum6)
mummodel6cri <- HPDinterval(mum6, 0.95)

dad6 <- model6$VCV[,"dad"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
dadmodel6 <- posterior.mode(dad6)
dadmodel6cri <- HPDinterval(dad6, 0.95)

res6 <- model6$VCV[,"units"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
resmodel6 <- posterior.mode(res6)
resmodel6cri <- HPDinterval(res6, 0.95)


#~~ repeatability

rep6 <- (model6$VCV[,"id"]+model6$VCV[,"animal"])/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+model6$VCV[,"U_PlateID"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
repmodel6 <- posterior.mode(rep6)
repmodel6cri <- HPDinterval(rep6, 0.95)


#~~ h2 without plate

addgen6woplate <- model6$VCV[,"animal"]/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
addgenmodel6woplate <- posterior.mode(addgen6woplate)
addgenmodel6criwoplate <- HPDinterval(addgen6woplate, 0.95)


#~~ repeatability without plate

rep6woplate <- (model6$VCV[,"id"]+model6$VCV[,"animal"])/
  (model6$VCV[,"animal"]+model6$VCV[,"id"]+
     model6$VCV[,"mum"]+model6$VCV[,"dad"]+model6$VCV[,"units"])
repmodel6woplate <- posterior.mode(rep6woplate)
repmodel6criwoplate <- HPDinterval(rep6woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Dad", "Residual")
VarEst6 <- rbind(addgenvar6, perenvvar6, platevar6, mumvar6, dadvar6, resvar6)
VarCri6 <- rbind(addgencri6, perenvcri6, platecri6, mumcri6, dadcri6, resvarcri6)
PropEst6 <- rbind(addgenmodel6, perenvmodel6, platemodel6, mummodel6, dadmodel6, resmodel6)
PropCri6 <- rbind(addgenmodel6cri, perenvmodel6cri, platemodel6cri, mummodel6cri, dadmodel6cri, resmodel6cri)

Results6 <- as.data.frame(cbind(Variance, VarEst6, VarCri6, PropEst6, PropCri6))
Results6$Model <- "6"


# # -----------------------------------------------------------------------

# Model 7

# # -----------------------------------------------------------------------


prior6 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G5=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G6=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model7<-MCMCglmm(RTL~Sex + logAge + Technician,
                 random=~animal + id + U_PlateID + mum + dad + FPID,
                 pedigree=ped,
                 data=telodata.ped,
                 nitt=3600000, thin=1500, burnin=600000, prior=prior6)

save(model7, file="AnMod_MCMCglmm_7_060320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_7_060320.RData")


#~~ diagnostics

plot(model7$Sol) 
plot(model7$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model7$Sol)
autocorr(model7$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model7$Sol)
effectiveSize(model7$VCV)


#~~ heidel

heidel.diag(model7$Sol)
heidel.diag(model7$VCV)


#~~ geweke 

geweke.diag(model7$VCV)
geweke.diag(model7$Sol)


#~~ variance estimates 

addgenvar7 <- posterior.mode(model7$VCV[,"animal"])
addgencri7 <- HPDinterval(model7$VCV[,"animal"], 0.95)

perenvvar7 <- posterior.mode(model7$VCV[,"id"])
perenvcri7 <- HPDinterval(model7$VCV[,"id"], 0.95)

platevar7 <- posterior.mode(model7$VCV[,"U_PlateID"])
platecri7 <- HPDinterval(model7$VCV[,"U_PlateID"], 0.95)

mumvar7 <- posterior.mode(model7$VCV[,"mum"])
mumcri7 <- HPDinterval(model7$VCV[,"mum"], 0.95)

dadvar7 <- posterior.mode(model7$VCV[,"dad"])
dadcri7 <- HPDinterval(model7$VCV[,"dad"], 0.95)

fpidvar7 <- posterior.mode(model7$VCV[,"FPID"])
fpidcri7 <- HPDinterval(model7$VCV[,"FPID"], 0.95)

resvar7 <- posterior.mode(model7$VCV[,"units"])
resvarcri7 <- HPDinterval(model7$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen7 <- model7$VCV[,"animal"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
addgenmodel7 <- posterior.mode(addgen7)
addgenmodel7cri <- HPDinterval(addgen7, 0.95)

perenv7 <- model7$VCV[,"id"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
perenvmodel7 <- posterior.mode(perenv7)
perenvmodel7cri <- HPDinterval(perenv7, 0.95)

plate7 <- model7$VCV[,"U_PlateID"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
platemodel7 <- posterior.mode(plate7)
platemodel7cri <- HPDinterval(plate7, 0.95)

mum7 <- model7$VCV[,"mum"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
mummodel7 <- posterior.mode(mum7)
mummodel7cri <- HPDinterval(mum7, 0.95)

dad7 <- model7$VCV[,"dad"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
dadmodel7 <- posterior.mode(dad7)
dadmodel7cri <- HPDinterval(dad7, 0.95)

fpid7 <- model7$VCV[,"FPID"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
fpidmodel7 <- posterior.mode(fpid7)
fpidmodel7cri <- HPDinterval(fpid7, 0.95)

res7 <- model7$VCV[,"units"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
resmodel7 <- posterior.mode(res7)
resmodel7cri <- HPDinterval(res7, 0.95)


#~~ repeatability

rep7 <- (model7$VCV[,"id"]+model7$VCV[,"animal"])/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"U_PlateID"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
repmodel7 <- posterior.mode(rep7)
repmodel7cri <- HPDinterval(rep7, 0.95)


#~~ h2 wo plate

addgen7woplate <- model7$VCV[,"animal"]/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
addgenmodel7woplate <- posterior.mode(addgen7woplate)
addgenmodel7criwoplate <- HPDinterval(addgen7woplate, 0.95)


#~~ rep wo plate

rep7woplate <- (model7$VCV[,"id"]+model7$VCV[,"animal"])/
  (model7$VCV[,"animal"]+model7$VCV[,"id"]+model7$VCV[,"FPID"]+
     model7$VCV[,"mum"]+model7$VCV[,"dad"]+model7$VCV[,"units"])
repmodel7woplate <- posterior.mode(rep7woplate)
repmodel7criwoplate <- HPDinterval(rep7woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Dad", "FPID", "Residual")
VarEst7 <- rbind(addgenvar7, perenvvar7, platevar7, mumvar7, dadvar7, fpidvar7, resvar7)
VarCri7 <- rbind(addgencri7, perenvcri7, platecri7, mumcri7, dadcri7, fpidcri7, resvarcri7)
PropEst7 <- rbind(addgenmodel7, perenvmodel7, platemodel7, mummodel7, dadmodel7, fpidmodel7, resmodel7)
PropCri7 <- rbind(addgenmodel7cri, perenvmodel7cri, platemodel7cri, mummodel7cri, dadmodel7cri, fpidmodel7cri, resmodel7cri)

Results7 <- as.data.frame(cbind(Variance, VarEst7, VarCri7, PropEst7, PropCri7))
Results7$Model <- "7"


# # -----------------------------------------------------------------------

# Model 8

# # -----------------------------------------------------------------------


prior7 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G5=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G6=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G7=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model8 <- MCMCglmm(RTL~Sex + logAge + Technician,
                   random=~animal + id + U_PlateID + mum + dad + FPID + Terr,
                   pedigree=ped,
                   data=telodata.ped,
                   nitt=3600000, thin=1500, burnin=600000, prior=prior7)

save(model8, file="AnMod_MCMCglmm_8_060320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_8_060320.RData")


#~~ diagnostics

plot(model8$Sol) 
plot(model8$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model8$Sol)
autocorr(model8$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model8$Sol)
effectiveSize(model8$VCV)


#~~ heidel

heidel.diag(model8$Sol)
heidel.diag(model8$VCV)


#~~ geweke 

geweke.diag(model8$VCV)
geweke.diag(model8$Sol)


#~~ variance estimates 

addgenvar8 <- posterior.mode(model8$VCV[,"animal"])
addgencri8 <- HPDinterval(model8$VCV[,"animal"], 0.95)

perenvvar8 <- posterior.mode(model8$VCV[,"id"])
perenvcri8 <- HPDinterval(model8$VCV[,"id"], 0.95)

platevar8 <- posterior.mode(model8$VCV[,"U_PlateID"])
platecri8 <- HPDinterval(model8$VCV[,"U_PlateID"], 0.95)

mumvar8 <- posterior.mode(model8$VCV[,"mum"])
mumcri8 <- HPDinterval(model8$VCV[,"mum"], 0.95)

dadvar8 <- posterior.mode(model8$VCV[,"dad"])
dadcri8 <- HPDinterval(model8$VCV[,"dad"], 0.95)

fpidvar8 <- posterior.mode(model8$VCV[,"FPID"])
fpidcri8 <- HPDinterval(model8$VCV[,"FPID"], 0.95)

terrvar8 <- posterior.mode(model8$VCV[,"Terr"])
terrcri8 <- HPDinterval(model8$VCV[,"Terr"], 0.95)

resvar8 <- posterior.mode(model8$VCV[,"units"])
resvarcri8 <- HPDinterval(model8$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen8 <- model8$VCV[,"animal"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
addgenmodel8 <- posterior.mode(addgen8)
addgenmodel8cri <- HPDinterval(addgen8, 0.95)

perenv8 <- model8$VCV[,"id"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
perenvmodel8 <- posterior.mode(perenv8)
perenvmodel8cri <- HPDinterval(perenv8, 0.95)

plate8 <- model8$VCV[,"U_PlateID"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
platemodel8 <- posterior.mode(plate8)
platemodel8cri <- HPDinterval(plate8, 0.95)

mum8 <- model8$VCV[,"mum"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
mummodel8 <- posterior.mode(mum8)
mummodel8cri <- HPDinterval(mum8, 0.95)

dad8 <- model8$VCV[,"dad"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
dadmodel8 <- posterior.mode(dad8)
dadmodel8cri <- HPDinterval(dad8, 0.95)

fpid8 <- model8$VCV[,"FPID"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
fpidmodel8 <- posterior.mode(fpid8)
fpidmodel8cri <- HPDinterval(fpid8, 0.95)

terr8 <- model8$VCV[,"Terr"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
terrmodel8 <- posterior.mode(terr8)
terrmodel8cri <- HPDinterval(terr8, 0.95)

res8 <- model8$VCV[,"units"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
resmodel8 <- posterior.mode(res8)
resmodel8cri <- HPDinterval(res8, 0.95)


#~~ repeatability

rep8 <- (model8$VCV[,"id"]+model8$VCV[,"animal"])/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"U_PlateID"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
repmodel8 <- posterior.mode(rep8)
repmodel8cri <- HPDinterval(rep8, 0.95)


#~~ h2 without plate

addgen8woplate <- model8$VCV[,"animal"]/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
addgenmodel8woplate <- posterior.mode(addgen8woplate)
addgenmodel8criwoplate <- HPDinterval(addgen8woplate, 0.95)


#~~ rep without plate 

rep8woplate <- (model8$VCV[,"id"]+model8$VCV[,"animal"])/
  (model8$VCV[,"animal"]+model8$VCV[,"id"]+model8$VCV[,"FPID"]+
     model8$VCV[,"mum"]+model8$VCV[,"dad"]+model8$VCV[,"units"]+model8$VCV[,"Terr"])
repmodel8woplate <- posterior.mode(rep8woplate)
repmodel8criwoplate <- HPDinterval(rep8woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Dad", "FPID", "Terr", "Residual")
VarEst8 <- rbind(addgenvar8, perenvvar8, platevar8, mumvar8, dadvar8, fpidvar8, terrvar8, resvar8)
VarCri8 <- rbind(addgencri8, perenvcri8, platecri8, mumcri8, dadcri8, fpidcri8, terrcri8, resvarcri8)
PropEst8 <- rbind(addgenmodel8, perenvmodel8, platemodel8, mummodel8, dadmodel8, fpidmodel8, terrmodel8, resmodel8)
PropCri8 <- rbind(addgenmodel8cri, perenvmodel8cri, platemodel8cri, mummodel8cri, dadmodel8cri, fpidmodel8cri, terrmodel8cri, resmodel8cri)

Results8 <- as.data.frame(cbind(Variance, VarEst8, VarCri8, PropEst8, PropCri8))
Results8$Model <- "8"


# # -----------------------------------------------------------------------

# Model 9

# # -----------------------------------------------------------------------


prior8 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G5=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G6=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G7=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G8=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model9<-MCMCglmm(RTL~Sex + logAge + Technician,
                 random=~animal + id + U_PlateID + mum + dad + FPID + Terr + BirthFPID,
                 pedigree=ped,
                 data=telodata.ped,
                 nitt=3600000, thin=1500, burnin=600000, prior=prior8)

save(model9, file="AnMod_MCMCglmm_9_060320.RData")

rm(list=ls()[! ls() %in% KeepThese])

load(file="AnMod_MCMCglmm_9_060320.RData")


#~~ diagnostics

plot(model9$Sol) 
plot(model9$VCV)


#~~ correlation among thinned samples must be low (<0.1)

autocorr(model9$Sol)
autocorr(model9$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model9$Sol)
effectiveSize(model9$VCV)


#~~ heidel

heidel.diag(model9$Sol)
heidel.diag(model9$VCV)


#~~ geweke 

geweke.diag(model9$VCV)
geweke.diag(model9$Sol)


#~~ fixed effects

posterior.mode(model9$Sol)
HPDinterval(model9$Sol, 0.95)


#~~ variance estimates 

addgenvar9 <- posterior.mode(model9$VCV[,"animal"])
addgencri9 <- HPDinterval(model9$VCV[,"animal"], 0.95)

perenvvar9 <- posterior.mode(model9$VCV[,"id"])
perenvcri9 <- HPDinterval(model9$VCV[,"id"], 0.95)

platevar9 <- posterior.mode(model9$VCV[,"U_PlateID"])
platecri9 <- HPDinterval(model9$VCV[,"U_PlateID"], 0.95)

mumvar9 <- posterior.mode(model9$VCV[,"mum"])
mumcri9 <- HPDinterval(model9$VCV[,"mum"], 0.95)

dadvar9 <- posterior.mode(model9$VCV[,"dad"])
dadcri9 <- HPDinterval(model9$VCV[,"dad"], 0.95)

fpidvar9 <- posterior.mode(model9$VCV[,"FPID"])
fpidcri9 <- HPDinterval(model9$VCV[,"FPID"], 0.95)

terrvar9 <- posterior.mode(model9$VCV[,"Terr"])
terrcri9 <- HPDinterval(model9$VCV[,"Terr"], 0.95)

bfpidvar9 <- posterior.mode(model9$VCV[,"BirthFPID"])
bfpidcri9 <- HPDinterval(model9$VCV[,"BirthFPID"], 0.95)

resvar9 <- posterior.mode(model9$VCV[,"units"])
resvarcri9 <- HPDinterval(model9$VCV[,"units"], 0.95)


#~~ evolvability

addgenvar9/(mean(telodata.ped$RTL)*mean(telodata.ped$RTL))
addgencri9/(mean(telodata.ped$RTL)*mean(telodata.ped$RTL))


#~~ prop of phenotypic variance explained

addgen9 <- model9$VCV[,"animal"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
addgenmodel9 <- posterior.mode(addgen9)
addgenmodel9cri <- HPDinterval(addgen9, 0.95)

perenv9 <- model9$VCV[,"id"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
perenvmodel9 <- posterior.mode(perenv9)
perenvmodel9cri <- HPDinterval(perenv9, 0.95)

plate9 <- model9$VCV[,"U_PlateID"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
platemodel9 <- posterior.mode(plate9)
platemodel9cri <- HPDinterval(plate9, 0.95)

mum9 <- model9$VCV[,"mum"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
mummodel9 <- posterior.mode(mum9)
mummodel9cri <- HPDinterval(mum9, 0.95)

dad9 <- model9$VCV[,"dad"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
dadmodel9 <- posterior.mode(dad9)
dadmodel9cri <- HPDinterval(dad9, 0.95)

fpid9 <- model9$VCV[,"FPID"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
fpidmodel9 <- posterior.mode(fpid9)
fpidmodel9cri <- HPDinterval(fpid9, 0.95)

terr9 <- model9$VCV[,"Terr"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
terrmodel9 <- posterior.mode(terr9)
terrmodel9cri <- HPDinterval(terr9, 0.95)

bfpid9 <- model9$VCV[,"BirthFPID"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
bfpidmodel9 <- posterior.mode(bfpid9)
bfpidmodel9cri <- HPDinterval(bfpid9, 0.95)

res9 <- model9$VCV[,"units"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
resmodel9 <- posterior.mode(res9)
resmodel9cri <- HPDinterval(res9, 0.95)


#~~ repeatability

rep9 <- (model9$VCV[,"id"]+model9$VCV[,"animal"])/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"U_PlateID"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
repmodel9 <- posterior.mode(rep9)
repmodel9cri <- HPDinterval(rep9, 0.95)

repmodel9
repmodel9cri


#~~ h2 wo plate

addgen9woplate <- model9$VCV[,"animal"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
addgenmodel9woplate <- posterior.mode(addgen9woplate)
addgenmodel9criwoplate <- HPDinterval(addgen9woplate, 0.95)

addgenmodel9woplate
addgenmodel9criwoplate


#~~ fpid wo plate

fpid9woplate <- model9$VCV[,"FPID"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
fpidmodel9woplate <- posterior.mode(fpid9woplate)
fpidmodel9criwoplate <- HPDinterval(fpid9woplate, 0.95)

fpidmodel9woplate
fpidmodel9criwoplate


#~~ repeatability wo plate

rep9woplate <- (model9$VCV[,"id"]+model9$VCV[,"animal"])/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
repmodel9woplate <- posterior.mode(rep9woplate)
repmodel9criwoplate <- HPDinterval(rep9woplate, 0.95)

repmodel9woplate
repmodel9criwoplate


#~~ res wo plate

res9woplate <- model9$VCV[,"units"]/
  (model9$VCV[,"animal"]+model9$VCV[,"id"]+model9$VCV[,"FPID"]+
     model9$VCV[,"mum"]+model9$VCV[,"dad"]+model9$VCV[,"units"]+model9$VCV[,"Terr"]+model9$VCV[,"BirthFPID"])
resmodel9woplate <- posterior.mode(res9woplate)
resmodel9criwoplate <- HPDinterval(res9woplate, 0.95)

resmodel9woplate
resmodel9criwoplate


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Dad", "FPID", "Terr", "BirthFPID", "Residual")
VarEst9 <- rbind(addgenvar9, perenvvar9, platevar9, mumvar9, dadvar9, fpidvar9, terrvar9, bfpidvar9, resvar9)
VarCri9 <- rbind(addgencri9, perenvcri9, platecri9, mumcri9, dadcri9, fpidcri9, terrcri9, bfpidcri9, resvarcri9)
PropEst9 <- rbind(addgenmodel9, perenvmodel9, platemodel9, mummodel9, dadmodel9, fpidmodel9, terrmodel9, bfpidmodel9, resmodel9)
PropCri9 <- rbind(addgenmodel9cri, perenvmodel9cri, platemodel9cri, mummodel9cri, dadmodel9cri, fpidmodel9cri, terrmodel9cri, bfpidmodel9cri, resmodel9cri)

Results9 <- as.data.frame(cbind(Variance, VarEst9, VarCri9, PropEst9, PropCri9))
Results9$Model <- "9"

FullResults19 <- rbind(Results, Results2, Results3,
                      Results4, Results5, Results6,
                     Results7, Results8, Results9)

FullResults19 <- FullResults19[,c(8,1:7)]

write.csv(FullResults19, file="AnModResults1to9_310320.csv", row.names = F)

rm(FullResults19)


# # -----------------------------------------------------------------------

# Graph for manuscript - figure 2 - animal model results 

# # -----------------------------------------------------------------------


FinalResults <- read.csv(file="AnModResults1to9_310320.csv", header = T)

varplot <- ggplot(FinalResults, aes(x = Model, y = var1.1, fill=factor(Variance, levels=c("Residual", "Plate", "BirthFPID", "Terr", "FPID", "Dad", "Mum", "PerEnv", "AddGen", "ID")))) + 
  geom_bar(stat='identity') +
  theme_classic() +
  scale_x_continuous(breaks = round(seq(1, 9, by = 1),1))  +
  scale_fill_manual(values=c("#fdb863","#e66101","#b2abd2","#5e3c99"), "Variance \ncomponent", 
                    limits=c("ID", "AddGen", "FPID", "Plate"), 
                    labels=c(bquote(V[ID]), bquote(V[A]), expression(V[CP]), expression(V[Plate]))) +
  ylab("Proportion of phenotypic variance") + xlab ("Model") + 
  theme(axis.title.x=element_text(colour="grey6", size=13),
        axis.title.y=element_text(colour="grey6", size=13),
        axis.line=element_line(colour="grey6"),
        text = element_text(size=12, colour="grey10"),
        legend.text=element_text(size=14, colour="grey6"),
        legend.title=element_text(colour="grey6", size=13),
        legend.text.align = 0) + 
  scale_y_continuous(expand = c(0, 0.01), limits = c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) 

varplot

pdf("Figure2.pdf", width=7, height=4)
print(varplot)
dev.off()


# # -----------------------------------------------------------------------

# Model 10 with dominant female and dominant male in natal territory included rather than genetic parents
# Same model structure as Model 7

# # -----------------------------------------------------------------------


rm(list=ls()[! ls() %in% KeepThese])

prior6 <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G5=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                      G6=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)),
               R=list(V=1, n=0.002))

model10<-MCMCglmm(RTL~Sex + logAge + Technician,
                  random=~animal + id + U_PlateID + BrF + BrM + FPID,
                  pedigree=ped,
                  data=telodata.ped,
                  nitt=10000000, thin=3000, burnin=4000000, prior=prior6)

save(model10, file="AnMod_MCMCglmm_10_BrFBrM_300620.RData")

load(file="AnMod_MCMCglmm_10_BrFBrM_300620.RData")


#~~ diagnostics

plot(model10$Sol) 
plot(model10$VCV)


#~~ correlation among thinned samples must be low (<0.1) 

autocorr(model10$Sol)
autocorr(model10$VCV)


#~~ also need to check effective sample sample - >1000 needed

effectiveSize(model10$Sol)
effectiveSize(model10$VCV)


#~~ heidel

heidel.diag(model10$Sol)
heidel.diag(model10$VCV)


#~~ geweke 

geweke.diag(model10$VCV)
geweke.diag(model10$Sol)


#~~ variance estimates 

addgenvar10 <- posterior.mode(model10$VCV[,"animal"])
addgencri10 <- HPDinterval(model10$VCV[,"animal"], 0.95)

perenvvar10 <- posterior.mode(model10$VCV[,"id"])
perenvcri10 <- HPDinterval(model10$VCV[,"id"], 0.95)

platevar10 <- posterior.mode(model10$VCV[,"U_PlateID"])
platecri10 <- HPDinterval(model10$VCV[,"U_PlateID"], 0.95)

mumvar10 <- posterior.mode(model10$VCV[,"BrF"])
mumcri10 <- HPDinterval(model10$VCV[,"BrF"], 0.95)

dadvar10 <- posterior.mode(model10$VCV[,"BrM"])
dadcri10 <- HPDinterval(model10$VCV[,"BrM"], 0.95)

fpidvar10 <- posterior.mode(model10$VCV[,"FPID"])
fpidcri10 <- HPDinterval(model10$VCV[,"FPID"], 0.95)

resvar10 <- posterior.mode(model10$VCV[,"units"])
resvarcri10 <- HPDinterval(model10$VCV[,"units"], 0.95)


#~~ prop of phenotypic variance explained

addgen10 <- model10$VCV[,"animal"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
addgenmodel10 <- posterior.mode(addgen10)
addgenmodel10cri <- HPDinterval(addgen10, 0.95)

perenv10 <- model10$VCV[,"id"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
perenvmodel10 <- posterior.mode(perenv10)
perenvmodel10cri <- HPDinterval(perenv10, 0.95)

plate10 <- model10$VCV[,"U_PlateID"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
platemodel10 <- posterior.mode(plate10)
platemodel10cri <- HPDinterval(plate10, 0.95)

mum10 <- model10$VCV[,"BrF"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
mummodel10 <- posterior.mode(mum10)
mummodel10cri <- HPDinterval(mum10, 0.95)

dad10 <- model10$VCV[,"BrM"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
dadmodel10 <- posterior.mode(dad10)
dadmodel10cri <- HPDinterval(dad10, 0.95)

fpid10 <- model10$VCV[,"FPID"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
fpidmodel10 <- posterior.mode(fpid10)
fpidmodel10cri <- HPDinterval(fpid10, 0.95)

res10 <- model10$VCV[,"units"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
resmodel10 <- posterior.mode(res10)
resmodel10cri <- HPDinterval(res10, 0.95)


#~~ repeatability

rep10 <- (model10$VCV[,"id"]+model10$VCV[,"animal"])/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"U_PlateID"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
repmodel10 <- posterior.mode(rep10)
repmodel10cri <- HPDinterval(rep10, 0.95)


#~~ h2 wo plate

addgen10woplate <- model10$VCV[,"animal"]/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
addgenmodel10woplate <- posterior.mode(addgen10woplate)
addgenmodel10criwoplate <- HPDinterval(addgen10woplate, 0.95)


#~~ rep wo plate

rep10woplate <- (model10$VCV[,"id"]+model10$VCV[,"animal"])/
  (model10$VCV[,"animal"]+model10$VCV[,"id"]+model10$VCV[,"FPID"]+
     model10$VCV[,"BrF"]+model10$VCV[,"BrM"]+model10$VCV[,"units"])
repmodel10woplate <- posterior.mode(rep10woplate)
repmodel10criwoplate <- HPDinterval(rep10woplate, 0.95)


#~~ for results

Variance <- c("AddGen", "PerEnv", "Plate", "Mum", "Dad", "FPID", "Residual")
VarEst10 <- rbind(addgenvar10, perenvvar10, platevar10, mumvar10, dadvar10, fpidvar10, resvar10)
VarCri10 <- rbind(addgencri10, perenvcri10, platecri10, mumcri10, dadcri10, fpidcri10, resvarcri10)
PropEst10 <- rbind(addgenmodel10, perenvmodel10, platemodel10, mummodel10, dadmodel10, fpidmodel10, resmodel10)
PropCri10 <- rbind(addgenmodel10cri, perenvmodel10cri, platemodel10cri, mummodel10cri, dadmodel10cri, fpidmodel10cri, resmodel10cri)

Results10 <- as.data.frame(cbind(Variance, VarEst10, VarCri10, PropEst10, PropCri10))
Results10$Model <- "10"

write.csv(Results10, file="AnModResults10BrFBrM_300620.csv", row.names = F)

  
# # -----------------------------------------------------------------------

# ASReml results for Model 9

# # -----------------------------------------------------------------------


#~~ pin function

pin.function <- function(model){
  x <- summary(model)$varcomp
  totvar   <- sum(x[,2])
  x$Effect <- x$component/totvar
  x$SE <- NA
  
  
  object <- model
  pframe <- as.list(object$gammas) 
  
  denominator <- "1 "
  ab<-length(pframe)-1
  
  for(effname in names(pframe)[c(1:ab) ]){
    denominator <- paste(denominator, "+ `", effname, "` ", sep = "")
  }
  
  denominator <- paste("(", denominator, ")", sep = "")
  
  
  for(effname in names(pframe)[c(1:ab)]){
    transform <- eval(parse(text = paste("`placeholder` ~ `", effname, "`/", denominator, sep = "")))
    
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe) 
    X <- as.vector(attr(tvalue, "gradient")) 
    tname <- if (length(transform) == 3) {
      transform[[2]] 
    } else "" 
    V <- object$ai 
    n <- length(pframe) 
    i <- rep(1:n, 1:n) 
    if (length(V) != length(i)) 
      stop("vcov matrix incompatible with\nnumber of variance components") 
    j <- sequence(1:n) 
    k <- 1 + (i > j) 
    se <- sqrt(sum(V * X[i] * X[j] * k))
    x[effname,"SE"] <- se
  }
  
  return(x)
}

ainv<-asreml.Ainverse(ped)$ginv


#~~ Model 9 in asreml

asreml.model9<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + U_PlateID + mum + dad + FPID + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

summary(asreml.model9)

save(asreml.model9, file="AnMod_ASReml_9_260320")

pin.function(asreml.model9)

varcompallages <- as.data.frame(pin.function(asreml.model9))

write.csv(varcompallages, file="ASReml_PropsWPlate_030420.csv", row.names = T)


#~~ SE for Resid - EST 0.671484 SE 0.03091493

Resid2<-nadiv:::pin(asreml.model9,resid~V9/(V1+V2+V3+V4+V5+V6+V7+V8+V9))

Resid2


#~~ repeatability - rep 0.05738527 SE 0.0234785

Rep2<-nadiv:::pin(asreml.model9,rep~(V7+V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9))

Rep2


#~~ now without plate

FPID<-nadiv:::pin(asreml.model9,fpid~V1/(V1+V2+V3+V5+V6+V7+V8+V9))
BFPID<-nadiv:::pin(asreml.model9,bfpid~V2/(V1+V2+V3+V5+V6+V7+V8+V9))
Terr<-nadiv:::pin(asreml.model9,terr~V3/(V1+V2+V3+V5+V6+V7+V8+V9))
Dad<-nadiv:::pin(asreml.model9,dad~V5/(V1+V2+V3+V5+V6+V7+V8+V9))
Mum<-nadiv:::pin(asreml.model9,mum~V6/(V1+V2+V3+V5+V6+V7+V8+V9))
AddGen<-nadiv:::pin(asreml.model9,addgen~V7/(V1+V2+V3+V5+V6+V7+V8+V9))
PerEnv<-nadiv:::pin(asreml.model9,perenv~V8/(V1+V2+V3+V5+V6+V7+V8+V9))
Resid<-nadiv:::pin(asreml.model9,resid~V9/(V1+V2+V3+V5+V6+V7+V8+V9))

AddGen


#~~ repeatability without plate

Rep<-nadiv:::pin(asreml.model9,rep~(V7+V8)/(V1+V2+V3+V5+V6+V7+V8+V9))
Rep

#~~ Estimate         SE
#~~ rep 0.07376325 0.03004592

FPID$list<-"FPID" 
BFPID$list<-"BFPID"
Terr$list<-"Terr"
Dad$list<-"Dad"
Mum$list<-"Mum"
AddGen$list<-"AddGen"
PerEnv$list<-"PerEnv"
Resid$list<-"Resid"

proportions<-rbind(FPID, BFPID, Terr, Dad, Mum, AddGen, PerEnv, Resid)

write.csv(proportions, file="ASReml_PropsWoPlate_030420.csv", row.names = F)


#~~ fixed effects

summary(asreml.model9, all=TRUE)$coef.fixed

wald.asreml(asreml.model9, ssType="conditional", denDF="numeric")


#~~ significance of random effects

#~~ add gen

asreml.model9a<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ide(id,var=T,init=1) + U_PlateID + mum + dad + FPID + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")



2*(asreml.model9$loglik-asreml.model9a$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9a$loglik),1)

#~~ addgen - x=7.914205, p=0.004904812


#~~ per env

asreml.model9b<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ U_PlateID + mum + dad + FPID + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")


2*(asreml.model9$loglik-asreml.model9b$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9b$loglik),1)

#~~ per env - x=0.5203222, p=0.4707043


#~~ plate

asreml.model9c<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + mum + dad + FPID + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

2*(asreml.model9$loglik-asreml.model9c$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9c$loglik),1)

#~~ plate - x=223.8918, p=0


#~~ mum

asreml.model9d<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + U_PlateID + dad + FPID + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

2*(asreml.model9$loglik-asreml.model9d$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9d$loglik),1)

#~~ mum - x=-52.71716, p= 1


#~~ fpid

asreml.model9e<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + U_PlateID + mum + dad + Terr + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

2*(asreml.model9$loglik-asreml.model9e$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9e$loglik),1)

#~~ fpid - x=17.53602, p=2.819165e-05


#~~ birth fpid

asreml.model9f<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + U_PlateID + mum + dad + FPID + Terr,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

2*(asreml.model9$loglik-asreml.model9f$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9f$loglik),1)

#~~ birth fpid - x=0.9395412  , p=0.3323959


#~~ terr

asreml.model9g<-asreml(fixed = RTL ~ Sex + logAge + Technician,
                       random =~ped(id,var=T,init=1)+ide(id,var=T,init=1) + U_PlateID + mum + dad + FPID + BirthFPID,
                       data=telodata.ped, ginverse=list(id=ainv), 
                       na.method.X="omit", na.method.Y="omit")

2*(asreml.model9$loglik-asreml.model9g$loglik)
1-pchisq(2*(asreml.model9$loglik-asreml.model9g$loglik),1)

#~~ terr - x=0.5162797, p=0.472433
