# Start -------------------------------------------------------------------
rm(list = ls(all=TRUE))
source("./functions.R")
library(geomorph)
library(ggplot2)
library(patchwork)
library(geiger)
library(tidyr)
library(dplyr)
library(phytools)
library(OUwie)
library(ratematrix)
library(parallel)

# Load Data ---------------------------------------------------------------

# read internal data
InData <- read.csv("InternalDataAllNew.csv")
rownames(InData) <- InData$Species
InData$Habitat <- factor(InData$Habitat , levels=c("Aq", "SAq", "T"))
InData$Life <- factor(InData$Life , levels=c("pd", "bi", "dd"))
InData$HabitatLife <- factor(InData$HabitatLife , levels=c("Aq_pd", "Aq_bi", "SAq_bi","T_bi", "T_dd"))
HumerusSub <- InData[,c(1:12)]
FemurSub <- na.omit(InData[,c(1:7,13:17)])

# read 3D GMM Data
Humerus_land <- read_SlicerLand("Humerus_Land/")
hum.gpa <- gpagen(Humerus_land)

Femur_land <- read_SlicerLand("Femur_Land/")
fem.gpa <- gpagen(Femur_land)

# read tree
tree <- read.tree("Supplementary_File_S3_fulltree_trim.tree")
hum.tree <- drop.tip(tree, setdiff(tree$tip.label,HumerusSub$Species))
fem.tree <- drop.tip(tree, setdiff(tree$tip.label,FemurSub$Species))

# set colors
col_habitat <- setNames(c("#0a97b7","#fade7d","#e76062"),c("Aq","SAq","T"))
col_life <- setNames(c("purple3","orange","forestgreen"),c("pd","bi","dd"))
col_habitatlife <- setNames(c("#0a97b7","#98dce8","#fade7d","#fea360", "#e76062"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"))
shape_habitatlife <- setNames(c(24,21,21,21,22), c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"))

# set variables
hum_SMA <- setNames(as.numeric(HumerusSub$Hum_SMA),rownames(HumerusSub))
hum_COM <- setNames(as.numeric(HumerusSub$Hum_COM),rownames(HumerusSub))
fem_SMA <- setNames(as.numeric(FemurSub$Fem_SMA),rownames(FemurSub))
fem_COM <- setNames(as.numeric(FemurSub$Fem_COM),rownames(FemurSub))

# Anc State Rec -----------------------------------------------------------

corHMM.data <- InData[c(1,5,4)]
plot(tree, show.tip.label = FALSE)
tiplabels(pch = 16, col = col_habitat[corHMM.data[,2]], cex = 0.75)
tiplabels(pch = 16, col = col_life[corHMM.data[,3]], cex = 0.75, offset = 2.5)

# ## no dual transitions
RunModels(phy = tree, data = corHMM.data, dual = FALSE, name = "corModelSet.StewartWiens.Rsave")

# ## if we want to allow dual transitions
RunModels(phy = tree, data = corHMM.data, dual = TRUE, name = "corModelSet.DUAL.StewartWiens.Rsave")

load("corModelSet.StewartWiens.Rsave")
CorRes_1 <- CorRes_i
load("corModelSet.DUAL.StewartWiens.Rsave")
CorRes_2 <- CorRes_i
CorRes_i <- c(CorRes_1, CorRes_2)
save(CorRes_i, file = "corModelSet.Full.StewartWiens.Rsave")

results <- getResultsTable("corModelSet.Full.StewartWiens.Rsave")
best_model <- which(results[,"AICcWt"] == max(results[,"AICcWt"]))

## best fitting model as determined by wtAIC
load("corModelSet.Full.StewartWiens.Rsave")
MK_5state <- CorRes_i[[best_model]]
MK_5state; plotMKmodel(MK_5state)

## states:
#1: 0_0 Aq_bi
#2: 0_2 Aq_pd
#3: 1_0 SAq_bi
#4: 2_0 T_bi
#5: 2_1 T_dd

# run simmap with best fitting model paramters
simmap <- makeSimmap(tree=tree, data=MK_5state$data, model=MK_5state$solution, rate.cat= MK_5state$rate.cat, nSim=1000)
saveRDS(simmap, file = "corHMM.simmap.StewartWiens.RDS") # save simmap
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS") # load simmap

## summarize simmap
pd <- describe.simmap(simmap)

# plot the simmap for figure 1
plotTree(ladderize(tree, right = T),type="fan",part=0.91, fsize = 0.01)
obj<-axis(1,pos=-6,at=seq(25,175,by=25),cex.axis=0.1,labels=FALSE, lwd = 1.5)
numb <- c(sort(obj, decreasing = T)[-1],0)
text(obj,rep(-18,length(obj)),numb,cex=0.8,lwd = 1.5)
text(mean(obj),-27,"time (mya)",cex=0.8)
for(i in 1:(length(obj)-1)){
  a1<-atan(-2/obj[i])
  a2<-0.88*2*pi
  plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1,
                    col=make.transparent("darkgrey",0.5))
}
nodelabels(pie = pd$ace[,c(2,1,3,4,5)], piecol = col_habitatlife, cex = 0.4)
tiplabels(pch = shape_habitatlife[InData[tree$tip.label,"HabitatLife"]], bg = col_habitatlife[InData[tree$tip.label,"HabitatLife"]], cex = 1.2)

# subset simmap to femur data
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS")
for (i in 1:length(simmap)) {
  simmap[[i]] <- drop.tip(simmap[[i]],setdiff(simmap[[i]]$tip.label, fem.tree$tip.label))
}
saveRDS(simmap,"fem.corHMM.simmap.StewartWiens.RDS")


# Plot 3D GMM -------------------------------------------------------------

# run PCAs
hum.pca <- gm.prcomp(hum.gpa$coords, phy = hum.tree, align.to.phy = F)
fem.pca <- gm.prcomp(fem.gpa$coords, phy = fem.tree, align.to.phy = F)


# plot the 3D GMM PCAs for figure 2
par(mfrow = c(2,2),mar = c(4, 4, 1, 1))
plot(hum.pca,1,2, bg = col_habitatlife[as.factor(HumerusSub$HabitatLife)],  flip = c(1,2),
     pch = c(24,21,22)[as.factor(HumerusSub$Life)], phylo = T, cex = 1.5,
     phylo.par= list(node.labels = F, anc.states = F, tip.labels = F, tip.txt.cex = 0.5,
                     edge.color = "grey40"))
plot(fem.pca,1,2, bg = col_habitatlife[as.factor(FemurSub$HabitatLife)],  flip = c(1,2),
     pch = c(24,21,22)[as.factor(FemurSub$Life)], phylo = T, cex = 1.5,
     phylo.par= list(node.labels = F, anc.states = F, tip.labels = F, tip.txt.cex = 0.5, 
                     edge.color = "grey40"))
plot(hum.pca,1,3, bg = col_habitatlife[as.factor(HumerusSub$HabitatLife)], flip = c(1),
     pch = c(24,21,22)[as.factor(HumerusSub$Life)], phylo = T, cex = 1.5, 
     phylo.par= list(node.labels = F, anc.states = F, tip.labels = F, tip.txt.cex = 0.5,
                     edge.color = "grey40"))
plot(fem.pca,1,3, bg = col_habitatlife[as.factor(FemurSub$HabitatLife)], flip = c(1,2),
     pch = c(24,21,22)[as.factor(FemurSub$Life)], phylo = T, cex = 1.5,
     phylo.par= list(node.labels = F, anc.states = F, tip.labels = F, tip.txt.cex = 0.5, 
                     edge.color = "grey40"))
par(mfrow=c(1,1))

# Plot Internal Traits ----------------------------------------------------

#HUMERUS SMA
hic <- ggplot(HumerusSub,aes(y = Hum_SMA, x = HabitatLife)) + 
  geom_boxplot(outlier.alpha = 0, aes(fill = HabitatLife))+
  geom_point(size = 2.5,position = position_jitter(0.25), alpha = 1, aes(shape = HabitatLife, fill = HabitatLife))+
  scale_shape_manual(values = c(24,21,21,21,22))+
  ylim(1,2.73)+
  labs(y = "Stiffness", x = "")+
  scale_fill_manual(values = col_habitatlife)+ theme_classic() 

#FEMUR SMA
fic <- ggplot(FemurSub,aes(y = Fem_SMA, x = HabitatLife, fill = HabitatLife,)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(size = 2.5,position = position_jitter(0.25), alpha = 1, aes(shape = HabitatLife,fill = HabitatLife))+
  scale_shape_manual(values = c(24,21,21,21,22))+
  ylim(1,2.73)+
  labs(y = "Stiffness")+
  scale_fill_manual(values = col_habitatlife)+ theme_classic() 

# HUMERUS COM
hcc <- ggplot(HumerusSub,aes(y = Hum_COM, x = HabitatLife, fill = HabitatLife)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(size = 2.5,position = position_jitter(0.25), alpha = 1, aes(shape = HabitatLife,fill = HabitatLife))+
  scale_shape_manual(values = c(24,21,21,21,22))+
  ylim(0.6,1.001)+
  labs(y = "Density")+
  scale_fill_manual(values = col_habitatlife)+ theme_classic() 

# FEMUR COM
fcc <- ggplot(FemurSub,aes(y = Fem_COM, x = HabitatLife, fill = HabitatLife)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(size = 2.5,position = position_jitter(0.25), alpha = 1, aes(shape = HabitatLife,fill = HabitatLife))+
  scale_shape_manual(values = c(24,21,21,21,22))+
  ylim(0.6,1.001)+
  labs(y = "Density")+
  scale_fill_manual(values = col_habitatlife)+ theme_classic() 

hic + fic + hcc + fcc & theme(legend.position = "none")


# ANOVAs ------------------------------------------------------------------
## External Shape --------------------------------------------------------
### Humerus GMM -----------------------------------------------------------

# with sirens
lm5 <- procD.pgls(hum.gpa$coords~log(HumerusSub$SVL)*(HabitatLife),data = HumerusSub, phy = hum.tree)
summary(lm5)

# to compare means between groups
null <-  procD.pgls(hum.gpa$coords~ log(HumerusSub$SVL), data = HumerusSub,phy = hum.tree)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife))

# to compare slopes...
null <-  procD.pgls(hum.gpa$coords~ log(HumerusSub$SVL)+(HabitatLife), data = HumerusSub,phy = hum.tree)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife, covariate = log(HumerusSub$SVL)))

allom <- plotAllometry(lm5, HumerusSub$SVL, logsz = T, method = "PredLine", pch = 21, 
                       bg = col_habitatlife[as.factor(HumerusSub$HabitatLife)])


# without sirens
# to compare means between groups
lm5 <- procD.pgls(hum.gpa$coords[,,FemurSub$Species]~log(SVL)*(HabitatLife),data = HumerusSub[FemurSub$Species,], phy = fem.tree)
summary(lm5)
null <-  procD.pgls(hum.gpa$coords[,,FemurSub$Species]~ log(SVL), data = HumerusSub[FemurSub$Species,],phy = fem.tree)
p <- summary(pairwise(lm5, null, groups = HumerusSub[FemurSub$Species,"HabitatLife"]))

# to compare slopes...
null <-  procD.pgls(hum.gpa$coords[,,FemurSub$Species]~ log(SVL)+(HabitatLife), data = HumerusSub[FemurSub$Species,],phy = fem.tree)
p <- summary(pairwise(lm5, null, groups = HumerusSub[FemurSub$Species,"HabitatLife"], covariate = log(HumerusSub[FemurSub$Species,"SVL"])))

allom <- plotAllometry(lm5, HumerusSub[FemurSub$Species,"SVL"], logsz = T, method = "CAC", pch = 21, 
                       bg = col_habitatlife[as.factor(HumerusSub[FemurSub$Species,"HabitatLife"])])


### Femur GMM -----------------------------------------------------------

lm5 <- procD.pgls(fem.gpa$coords~log(FemurSub$SVL)*(HabitatLife),data = FemurSub, phy = fem.tree)

# to compare means between groups
summary(lm5)
null <-  procD.pgls(fem.gpa$coords~ log(FemurSub$SVL), data = FemurSub,phy = fem.tree)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife))

# to compare slopes...
null <-  procD.pgls(fem.gpa$coords~ log(FemurSub$SVL)+(HabitatLife), data = FemurSub,phy = fem.tree)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife, covariate = log(FemurSub$SVL)))

allom <- plotAllometry(lm5, FemurSub$SVL, logsz = T, method = "PredLine", pch = 21, 
                       bg = col_habitatlife[as.factor(FemurSub$HabitatLife)])

## Internal Traits --------------------------------------------------------
### Humerus SMA --------------------------------------------------------

# check fit of BM, OU, and EB
hum.tree.mod <- hum.tree
hum.tree.mod$edge.length <- hum.tree.mod$edge.length * 100
bm <- fitContinuous(hum.tree.mod, hum_SMA, model = "BM")
ou <- fitContinuous(hum.tree.mod, hum_SMA, model = "OU")
eb <- fitContinuous(hum.tree.mod, hum_SMA, model = "EB")
aic.vals<-setNames(c(bm$opt$aicc,ou$opt$aicc,eb$opt$aicc),c("BM","OU","EB")); aic.vals
aic.w(aic.vals) # OU best
OUvcv <- vcv(corMartins(1, phy = hum.tree.mod, form = ~hum.tree$tip.label)) 

# with sirens
lm5 <- procD.pgls(hum_SMA~log(SVL)*(HabitatLife),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
summary(lm5)

# to compare the means...
null <- procD.pgls(hum_SMA~log(SVL),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(hum_SMA~log(SVL)+(HabitatLife),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife, covariate = log(HumerusSub$SVL)));p

plotAllometry(lm5, HumerusSub$SVL, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(HumerusSub$HabitatLife)])

# without sirens
OUvcv <- vcv(corMartins(1, phy = fem.tree, form = ~fem.tree$tip.label))
lm5 <- procD.pgls(hum_SMA[FemurSub$Species]~log(SVL)*(HabitatLife),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
summary(lm5)

# to compare the means...
null <- procD.pgls(hum_SMA[FemurSub$Species]~log(SVL),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(hum_SMA[FemurSub$Species]~log(SVL)+(HabitatLife),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife, covariate = log(FemurSub$SVL)))

plotAllometry(lm5, FemurSub$SVL, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(FemurSub$HabitatLife)])

### Humerus COM -----------------------------------------------------------

# check fit of BM, OU, and EB
bm <- fitContinuous(hum.tree.mod, hum_COM, model = "BM")
ou <- fitContinuous(hum.tree.mod, hum_COM, model = "OU")
eb <- fitContinuous(hum.tree.mod, hum_COM, model = "EB")
aic.vals<-setNames(c(bm$opt$aicc,ou$opt$aicc,eb$opt$aicc),c("BM","OU","EB")); aic.vals
aic.w(aic.vals) # OU best
OUvcv <- vcv(corMartins(1, phy = hum.tree, form = ~hum.tree$tip.label))

# with sirens
lm5 <- procD.pgls(hum_COM~log(SVL)*(HabitatLife),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
summary(lm5)

# to compare the means
null <- procD.pgls(hum_COM~log(SVL),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(hum_COM~log(SVL)+(Habitat+Life),data = HumerusSub,phy = hum.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = HumerusSub$HabitatLife, covariate = log(HumerusSub$SVL)))

plotAllometry(lm5, HumerusSub$SVL, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(HumerusSub$HabitatLife)])

# without sirens
OUvcv <- vcv(corMartins(1, phy = fem.tree, form = ~fem.tree$tip.label))
lm5 <- procD.pgls(hum_COM[FemurSub$Species]~log(SVL)*(HabitatLife),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
summary(lm5)

# to compare the means...
null <- procD.pgls(hum_COM[FemurSub$Species]~log(SVL),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(hum_COM[FemurSub$Species]~log(SVL)+(HabitatLife),data = HumerusSub[FemurSub$Species,],phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife, covariate = log(FemurSub$SVL)))

plotAllometry(lm5, FemurSub$SVL, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(FemurSub$HabitatLife)])


### Femur SMA -----------------------------------------------------------

# check fit of BM, OU, and EB
fem.tree.mod <- fem.tree
fem.tree.mod$edge.length <- fem.tree.mod$edge.length * 100
bm <- fitContinuous(fem.tree.mod, fem_SMA, model = "BM")
ou <- fitContinuous(fem.tree.mod, fem_SMA, model = "OU")
eb <- fitContinuous(fem.tree.mod, fem_SMA, model = "EB")
aic.vals<-setNames(c(bm$opt$aicc,ou$opt$aicc,eb$opt$aicc),c("BM","OU","EB")); aic.vals
aic.w(aic.vals) # OU best
OUvcv <- vcv(corMartins(1, phy = fem.tree, form = ~fem.tree$tip.label))

lm5 <- procD.pgls(fem_SMA~log(SVL)*(HabitatLife),data = FemurSub,phy = fem.tree, Cov = OUvcv)
summary(lm5)

# to compare the means...
null <- procD.pgls(fem_SMA~log(SVL),data = FemurSub,phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(fem_SMA~log(SVL)+(Habitat+Life),data = FemurSub,phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife, covariate = log(FemurSub$SVL)))

plotAllometry(lm5, FemurSub$SVL, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(FemurSub$HabitatLife)])

### Femur COM -----------------------------------------------------------

# check fit of BM, OU, and EB
bm <- fitContinuous(fem.tree.mod, fem_COM, model = "BM")
ou <- fitContinuous(fem.tree.mod, fem_COM, model = "OU")
eb <- fitContinuous(fem.tree.mod, fem_COM, model = "EB")
aic.vals<-setNames(c(bm$opt$aicc,ou$opt$aicc,eb$opt$aicc),c("BM","OU","EB")); aic.vals
aic.w(aic.vals) # OU best
OUvcv <- vcv(corMartins(1, phy = fem.tree, form = ~fem.tree$tip.label))

lm5 <- procD.pgls(fem_COM~log(SVL)*(HabitatLife),data = FemurSub,phy = fem.tree, Cov = OUvcv)
summary(lm5)

# to compare the means...
null <- procD.pgls(fem_COM~log(SVL),data = FemurSub,phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife))

# to compare slopes...
null <- procD.pgls(fem_COM~log(SVL)+(HabitatLife),data = FemurSub,phy = fem.tree, Cov = OUvcv)
p <- summary(pairwise(lm5, null, groups = FemurSub$HabitatLife, covariate = log(FemurSub$SVL)))

plotAllometry(lm5, fem.gpa$Csize, logsz = T, method = "PredLine", pch = 21, 
              bg = col_habitatlife[as.factor(FemurSub$HabitatLife)])

## PairedTest -------------------------------------------------------------

# perform paired T test to compare mean humerus and femur trait by ecomorph

# SMA
ecotype <- levels(InData$HabitatLife)
for (e in 1:5) {
  group <- FemurSub$Species[which(FemurSub$HabitatLife == ecotype[e])]
  phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,group))
  
  hsma <- hum_SMA[phy$tip.label]
  hcom <- hum_COM[phy$tip.label]
  fsma <- fem_SMA[phy$tip.label]
  fcom <- fem_COM[phy$tip.label]
  
  print(ecotype[e])
  print(phyl.pairedttest(phy, hsma,fsma))
}

# COM
ecotype <- levels(InData$HabitatLife)
for (e in 1:5) {
  group <- FemurSub$Species[which(FemurSub$HabitatLife == ecotype[e])]
  phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,group))
  
  hsma <- hum_SMA[phy$tip.label]
  hcom <- hum_COM[phy$tip.label]
  fsma <- fem_SMA[phy$tip.label]
  fcom <- fem_COM[phy$tip.label]
  
  print(ecotype[e])
  print(phyl.pairedttest(phy, hcom,fcom))
}


# Disparity ---------------------------------------------------------------
## External Shape --------------------------------------------------------
### Humerus GMM --------------------------------------------------------

#with sirens
dis <- morphol.disparity(hum.gpa$coords~HumerusSub$HabitatLife, groups = HumerusSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

#without sirens
dis <- morphol.disparity(hum.gpa$coords[,,FemurSub$Species]~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

### Femur GMM --------------------------------------------------------

dis <- morphol.disparity(fem.gpa$coords~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/(dis$Procrustes.var)[3],3)

## Internal Traits --------------------------------------------------------
### Humerus SMA --------------------------------------------------------

# with sirens
dis <- morphol.disparity(hum_SMA~HumerusSub$HabitatLife, groups = HumerusSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

#without sirens
dis <- morphol.disparity(hum_SMA[FemurSub$Species]~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

### Humerus COM --------------------------------------------------------

# with sirens
dis <- morphol.disparity(hum_COM~HumerusSub$HabitatLife, groups = HumerusSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

# without sirens
dis <- morphol.disparity(hum_COM[FemurSub$Species]~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

### Femur SMA --------------------------------------------------------

dis <- morphol.disparity(fem_SMA~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

### Femur COM --------------------------------------------------------

dis <- morphol.disparity(fem_COM~FemurSub$HabitatLife, groups = FemurSub$HabitatLife);dis
dis$Procrustes.var
round(dis$Procrustes.var/min(dis$Procrustes.var),3)

# External Evo Rates -------------------------------------------------------

## Humerus GMM -------------------------------------------------------------

# with sirens
x <- compare.evol.rates(hum.gpa$coords, gp = setNames(HumerusSub$HabitatLife,HumerusSub$Species), phy = hum.tree)
round(x$sigma.d.gp,9)

# without sirens
x <- compare.evol.rates(hum.gpa$coords[,,FemurSub$Species], gp = setNames(FemurSub$HabitatLife,FemurSub$Species), phy = fem.tree)
round(x$sigma.d.gp,9)

## Femur GMM -------------------------------------------------------------

y <- compare.evol.rates(fem.gpa$coords, gp = setNames(FemurSub$HabitatLife,FemurSub$Species), phy = fem.tree)
round(y$sigma.d.gp,9)

# Ratematrix --------------------------------------------------------------
## Humerus -----------------------------------------------------------------
# testing if evolutionary correlations  differ between ecotypes
h.resp.data <- data.frame(row.names=names(hum_SMA),
                          hum_SMA=hum_SMA,
                          hum_COM=hum_COM)
h.pred.data <- setNames(HumerusSub[hum.tree$tip.label,"HabitatLife"],hum.tree$tip.label)

# explore distribution of predictor
summary( h.pred.data ) / length(h.pred.data)

# use the ASR and Q matrices to set priors
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS")
fit_Q_sym <- simmap[[1]]$Q
fit_Q_sym <- fit_Q_sym[c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd")]

maps <- lapply(1:10, function(x) fastSimmap(tree = hum.tree, x = h.pred.data, Q = fit_Q_sym))
## Now using a simple for loop.
maps <- vector(mode = "list", length = 10)
for( i in 1:1000 ) maps[[i]] <- fastSimmap(tree = hum.tree, x = h.pred.data, Q = fit_Q_sym)

# Run MCMC in parallel
# Define the number of cores to use
num_cores<-6
# Create a cluster object using the detected cores
cl <- makeCluster(num_cores)

clusterCall(cl, ratematrixMCMC, data = h.resp.data, 
            phy = maps, gen = 5000000, dir = "hum_runs_ASR_StewartWiens",
            burn = 0.25)
# Stop the cluster when done
stopCluster(cl)

# read in and check the results
files <- list.files("hum_runs_ASR_StewartWiens",full.names = T)
files <- files[grep(".rds",files)]
nruns <- length(files)
h.runlist <- list()
h.mcmc <- list()
h.corrpost <- list()
par(mfrow = c(ceiling(nruns/5),5))
for (i in 1:nruns) {
  h.runlist <- c(h.runlist,list(readRDS(files[i])))
  h.mcmc <- c(h.mcmc,list(readMCMC(h.runlist[[i]])))
  h.corrpost <- c(h.corrpost,list(extractCorrelation(post = h.mcmc[[i]])))
  computeESS(h.mcmc[[i]], p = 3)
  logAnalyzer(h.runlist[[i]])
  boxplot(h.corrpost[[i]], col = col_habitatlife, ylim = c(-1,1))
}

# check for convergence
checkConvergence(h.mcmc)

# save the merged posterior distrubitions
h.merged.mcmc <- mergePosterior(h.mcmc)
saveRDS(h.merged.mcmc,"hum_merged_mcmc.StewartWiens.RDS")


## Humerus no sirens -------------------------------------------------------
# testing if evolutionary correlations differ between ecotypes
h.resp.data <- data.frame(row.names=names(hum_SMA[fem.tree$tip.label]),
                          hum_SMA=hum_SMA[fem.tree$tip.label],
                          hum_COM=hum_COM[fem.tree$tip.label])
h.pred.data <- setNames(HumerusSub[fem.tree$tip.label,"HabitatLife"],fem.tree$tip.label)

# explore distribution of predictor
summary(h.pred.data) / length(h.pred.data)

# use the ASR and Q matrices to set priors
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS")
fit_Q_sym <- simmap[[1]]$Q
fit_Q_sym <- fit_Q_sym[c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd")]

maps <- lapply(1:10, function(x) fastSimmap(tree = fem.tree, x = h.pred.data, Q = fit_Q_sym))
## Now using a simple for loop.
maps <- vector(mode = "list", length = 10)
for( i in 1:1000 ) maps[[i]] <- fastSimmap(tree = fem.tree, x = h.pred.data, Q = fit_Q_sym)

# Run MCMC in parallel
# Define the number of cores to use
num_cores<-6
# Create a cluster object using the detected cores
cl <- makeCluster(num_cores)

clusterCall(cl, ratematrixMCMC, data = h.resp.data,
            phy = maps, gen = 5000000, dir = "hum_runs_nosirens_ASR_StewartWiens",
            burn = 0.25)
# Stop the cluster when done
stopCluster(cl)

# read in and check the results
files <- list.files("hum_runs_nosirens_ASR_StewartWiens",full.names = T)
files <- files[grep(".rds",files)]
nruns <- length(files)
h.runlist <- list()
h.mcmc <- list()
h.corrpost <- list()
par(mfrow = c(ceiling(nruns/5),5))
for (i in 1:nruns) {
  h.runlist <- c(h.runlist,list(readRDS(files[i])))
  h.mcmc <- c(h.mcmc,list(readMCMC(h.runlist[[i]])))
  h.corrpost <- c(h.corrpost,list(extractCorrelation(post = h.mcmc[[i]])))
  computeESS(h.mcmc[[i]], p = 3)
  logAnalyzer(h.runlist[[i]])
  boxplot(h.corrpost[[i]], col = col_habitatlife, ylim = c(-1,1))
}

# check for convergence
checkConvergence(h.mcmc)

# save the merged posterior distrubitions
h.merged.mcmc <- mergePosterior(h.mcmc)
saveRDS(h.merged.mcmc,"hum_nosirens_merged_mcmc.StewartWiens.RDS")


## Femur -------------------------------------------------------------------
# testing if evolutionary correlations differ between regimes
f.resp.data <- data.frame(row.names=names(fem_SMA),
                          fem_SMA=fem_SMA,
                          fem_COM=fem_COM)
f.pred.data <- setNames(FemurSub[fem.tree$tip.label,"HabitatLife"],fem.tree$tip.label)

# explore distribution of predictor
summary(f.pred.data) / length(f.pred.data)

# use the ASR and Q matrices to set priors
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS")
fit_Q_sym <- simmap[[1]]$Q
fit_Q_sym <- fit_Q_sym[c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd")]

maps <- lapply(1:10, function(x) fastSimmap(tree = fem.tree, x = f.pred.data, Q = fit_Q_sym))
## Now using a simple for loop.
maps <- vector(mode = "list", length = 10)
for( i in 1:1000 ) maps[[i]] <- fastSimmap(tree = fem.tree, x = f.pred.data, Q = fit_Q_sym)

# Run MCMC in parallel
# Define the number of cores to use
num_cores<-6
# Create a cluster object using the detected cores
cl <- makeCluster(num_cores)

clusterCall(cl, ratematrixMCMC, data = f.resp.data,
            phy = maps, gen = 5000000, dir = "fem_runs_ASR_StewartWiens",
            burn = 0.25)
# Stop the cluster when done
stopCluster(cl)

# read in and check the results
files <- list.files("fem_runs_ASR_StewartWiens",full.names = T)
files <- files[grep(".rds",files)]
nruns <- length(files)
f.runlist <- list()
f.mcmc <- list()
f.corrpost <- list()
par(mfrow = c(ceiling(nruns/5),5))
for (i in 1:nruns) {
  f.runlist <- c(f.runlist,list(readRDS(files[i])))
  f.mcmc <- c(f.mcmc,list(readMCMC(f.runlist[[i]])))
  f.corrpost <- c(f.corrpost,list(extractCorrelation(post = f.mcmc[[i]])))
  computeESS(h.mcmc[[i]], p = 3)
  logAnalyzer(h.runlist[[i]])
  boxplot(f.corrpost[[i]], col = col_habitatlife, ylim = c(-1,1))
}

# check for convergence
checkConvergence(f.mcmc)

# save the merged posterior distrubitions
f.merged.mcmc <- mergePosterior(f.mcmc)
saveRDS(f.merged.mcmc,"fem_merged_mcmc.StewartWiens.RDS")


## SMA --------------------------------------------------------------------
# testing if evolutionary correlations differ between regimes
sma.resp.data <- data.frame(row.names=names(fem_SMA[fem.tree$tip.label]),
                            hum_SMA=hum_SMA[fem.tree$tip.label],
                            fem_SMA=fem_SMA[fem.tree$tip.label])
sma.pred.data <- setNames(FemurSub[fem.tree$tip.label,"HabitatLife"],fem.tree$tip.label)

# explore distribution of predictor
summary(sma.pred.data) / length(sma.pred.data)

# use the ASR and Q matrices to set priors
simmap <- readRDS("corHMM.simmap.JetzPyron.RDS")
fit_Q_sym <- simmap[[1]]$Q
fit_Q_sym <- fit_Q_sym[c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd")]

maps <- lapply(1:10, function(x) fastSimmap(tree = fem.tree, x = sma.pred.data, Q = fit_Q_sym))
## Now using a simple for loop.
maps <- vector(mode = "list", length = 10)
for( i in 1:1000 ) maps[[i]] <- fastSimmap(tree = fem.tree, x = sma.pred.data, Q = fit_Q_sym)

# Run MCMC in parallel
# Define the number of cores to use
num_cores<-6
# Create a cluster object using the detected cores
cl <- makeCluster(num_cores)

clusterCall(cl, ratematrixMCMC, data = sma.resp.data,
            phy = maps, gen = 5000000, dir = "sma_runs_ASR_StewartWiens",
            burn = 0.25)
# Stop the cluster when done
stopCluster(cl)

# read in and check the results
files <- list.files("sma_runs_ASR_StewartWiens",full.names = T)
files <- files[grep(".rds",files)]
nruns <- length(files)
sma.runlist <- list()
sma.mcmc <- list()
sma.corrpost <- list()
par(mfrow = c(ceiling(nruns/5),5))
for (i in 1:nruns) {
  sma.runlist <- c(sma.runlist,list(readRDS(files[i])))
  sma.mcmc <- c(sma.mcmc,list(readMCMC(sma.runlist[[i]])))
  sma.corrpost <- c(sma.corrpost,list(extractCorrelation(post = sma.mcmc[[i]])))
  computeESS(h.mcmc[[i]], p = 3)
  logAnalyzer(h.runlist[[i]])
  boxplot(sma.corrpost[[i]], col = col_habitatlife, ylim = c(-1,1))
}

# check for convergence
checkConvergence(sma.mcmc)

# save the merged posterior distrubitions
sma.merged.mcmc <- mergePosterior(sma.mcmc)
saveRDS(sma.merged.mcmc,"sma_merged_mcmc.StewartWiens.RDS")

## COM --------------------------------------------------------------------
# testing if evolutionary correlations differ between regimes
com.resp.data <- data.frame(row.names=names(fem_COM[fem.tree$tip.label]),
                            hum_COM=hum_COM[fem.tree$tip.label],
                            fem_COM=fem_COM[fem.tree$tip.label])
com.pred.data <- setNames(FemurSub[fem.tree$tip.label,"HabitatLife"],fem.tree$tip.label)

# explore distribution of predictor
summary(com.pred.data) / length(com.pred.data)

# use the ASR and Q matrices to set priors
simmap <- readRDS("corHMM.simmap.StewartWiens.RDS")
fit_Q_sym <- simmap[[1]]$Q
fit_Q_sym <- fit_Q_sym[c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd"),c("Aq_pd","Aq_bi","SAq_bi","T_bi","T_dd")]

maps <- lapply(1:10, function(x) fastSimmap(tree = fem.tree, x = com.pred.data, Q = fit_Q_sym))
## Now using a simple for loop.
maps <- vector(mode = "list", length = 10)
for( i in 1:1000 ) maps[[i]] <- fastSimmap(tree = fem.tree, x = com.pred.data, Q = fit_Q_sym)

# Run MCMC in parallel
# Define the number of cores to use
num_cores<-6
# Create a cluster object using the detected cores
cl <- makeCluster(num_cores)

clusterCall(cl, ratematrixMCMC, data = com.resp.data,
            phy = maps, gen = 5000000, dir = "com_runs_ASR_StewartWiens",
            burn = 0.25)
# Stop the cluster when done
stopCluster(cl)

# read in and check the results
files <- list.files("com_runs_ASR_StewartWiens",full.names = T)
files <- files[grep(".rds",files)]
nruns <- length(files)
com.runlist <- list()
com.mcmc <- list()
com.corrpost <- list()
par(mfrow = c(ceiling(nruns/5),5))
for (i in 1:nruns) {
  com.runlist <- c(com.runlist,list(readRDS(files[i])))
  com.mcmc <- c(com.mcmc,list(readMCMC(com.runlist[[i]])))
  com.corrpost <- c(com.corrpost,list(extractCorrelation(post = com.mcmc[[i]])))
  computeESS(h.mcmc[[i]], p = 3)
  logAnalyzer(h.runlist[[i]])
  boxplot(com.corrpost[[i]], col = col_habitatlife, ylim = c(-1,1))
}

# check for convergence
checkConvergence(com.mcmc[1:5])

# save the merged posterior distrubitions
com.merged.mcmc <- mergePosterior(com.mcmc)
saveRDS(com.merged.mcmc,"com_merged_mcmc.StewartWiens.RDS")

## Plot Correlations ------------------------------------------------------

h.merged.mcmc <- readRDS("hum_merged_mcmc.StewartWiens.RDS")
h.merged.mcmc.nosirens <- readRDS("hum_nosirens_merged_mcmc.StewartWiens.RDS")
f.merged.mcmc <- readRDS("fem_merged_mcmc.StewartWiens.RDS")
sma.merged.mcmc <- readRDS("sma_merged_mcmc.StewartWiens.RDS")
com.merged.mcmc <- readRDS("com_merged_mcmc.StewartWiens.RDS")

h.corr <- extractCorrelation(h.merged.mcmc)
h.corr <- extractCorrelation(h.merged.mcmc.nosirens)
f.corr <- extractCorrelation(f.merged.mcmc)
sma.corr <- extractCorrelation(sma.merged.mcmc)
com.corr <- extractCorrelation(com.merged.mcmc)

# calculate mean correlation values
apply(h.corr,2,mean)
apply(h.corr.nosirens,2,mean)
apply(f.corr,2,mean)
apply(sma.corr,2,mean)
apply(com.corr,2,mean)

# plot the posterior distributions for figure 4 
corr <- data.frame("HabitatLife" = c(rep("Aq_pd",nrow(h.corr)),rep("Aq_bi",nrow(h.corr)),rep("SAq_bi",nrow(h.corr)),rep("T_bi",nrow(h.corr)),rep("T_dd",nrow(h.corr))),
                   "Correlation" = c(h.corr[,1],h.corr[,2],h.corr[,3],h.corr[,4],h.corr[,5]))
a<-ggplot(corr, aes(x=Correlation, fill=HabitatLife)) +
  geom_density(alpha=0.7, adjust = 2.5) +
  scale_fill_manual(values = col_habitatlife)+
  xlim(-1,1)+
  xlab(label = "Humerus: Stiffness vs Density")+
  theme_classic()

corr <- data.frame("HabitatLife" = c(rep("Aq_pd",nrow(f.corr)),rep("Aq_bi",nrow(f.corr)),rep("SAq_bi",nrow(f.corr)),rep("T_bi",nrow(f.corr)),rep("T_dd",nrow(f.corr))),
                   "Correlation" = c(f.corr[,1],f.corr[,2],f.corr[,3],f.corr[,4],f.corr[,5]))
b<-ggplot(corr, aes(x=Correlation, fill=HabitatLife)) +
  geom_density(alpha=0.7, adjust = 2.5) +
  scale_fill_manual(values = col_habitatlife)+
  xlim(-1,1)+
  xlab(label = "Femur: Stiffness vs Density")+
  theme_classic()

corr <- data.frame("HabitatLife" = c(rep("Aq_pd",nrow(sma.corr)),rep("Aq_bi",nrow(sma.corr)),rep("SAq_bi",nrow(sma.corr)),rep("T_bi",nrow(sma.corr)),rep("T_dd",nrow(sma.corr))),
                   "Correlation" = c(sma.corr[,1],sma.corr[,2],sma.corr[,3],sma.corr[,4],sma.corr[,5]))
c<-ggplot(corr, aes(x=Correlation, fill=HabitatLife)) +
  geom_density(alpha=0.7, adjust = 2.5) +
  scale_fill_manual(values = col_habitatlife)+
  xlim(-1,1)+
  xlab(label = "Stiffness: Humerus vs Femur")+
  theme_classic()

corr <- data.frame("HabitatLife" = c(rep("Aq_pd",nrow(com.corr)),rep("Aq_bi",nrow(com.corr)),rep("SAq_bi",nrow(com.corr)),rep("T_bi",nrow(com.corr)),rep("T_dd",nrow(com.corr))),
                   "Correlation" = c(com.corr[,1],com.corr[,2],com.corr[,3],com.corr[,4],com.corr[,5]))
d<-ggplot(corr, aes(x=Correlation, fill=HabitatLife)) +
  geom_density(alpha=0.7, adjust = 3) +
  scale_fill_manual(values = col_habitatlife)+
  xlim(-1,1)+
  xlab(label = "Density: Humerus vs Femur")+
  theme_classic()

(a+b) / (c+d) + plot_annotation(tag_levels = "A") & theme(legend.position = 'none') 

# Two Block PLS -----------------------------------------------------------

# integration of the external shapes between limbs
# external shape by ecotype
pd <- FemurSub$Species[which(FemurSub$HabitatLife == "Aq_pd")]
IT_pd <- phylo.integration(hum.gpa$coords[,,pd],
                           fem.gpa$coords[,,pd], phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,pd)))
summary(IT_pd)

aq_bi <- FemurSub$Species[which(FemurSub$HabitatLife == "Aq_bi")]
IT_aq_bi <- phylo.integration(hum.gpa$coords[,,aq_bi],
                           fem.gpa$coords[,,aq_bi], phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,aq_bi)))
summary(IT_aq_bi)

saq_bi <- FemurSub$Species[which(FemurSub$HabitatLife == "SAq_bi")]
IT_saq_bi <- phylo.integration(hum.gpa$coords[,,saq_bi],
                           fem.gpa$coords[,,saq_bi], phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,saq_bi)))
summary(IT_saq_bi)

t_bi <- FemurSub$Species[which(FemurSub$HabitatLife == "T_bi")]
IT_t_bi <- phylo.integration(hum.gpa$coords[,,t_bi],
                           fem.gpa$coords[,,t_bi], phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,t_bi)))
summary(IT_t_bi)

dd <- FemurSub$Species[which(FemurSub$HabitatLife == "T_dd")]
IT_dd <- phylo.integration(hum.gpa$coords[,,dd],
                           fem.gpa$coords[,,dd], phy = drop.tip(fem.tree,setdiff(fem.tree$tip.label,dd)))
summary(IT_dd)

# hOUwie -------------------------------------------------------------------

# create model structures
{models <- list(getOUParamStructure(model = "BM1", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OU1", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 2, null.model = T),               
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 1, null.model = F),
                getOUParamStructure(model = "BMV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUV", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUM", nObsState = 5, rate.cat = 2, null.model = T),
                getOUParamStructure(model = "OUMV", nObsState = 5, rate.cat = 2, null.model = T))
models <- setNames(models, c("cid_bm1","cid_ou1","cd_bmv","cd_ouv","cd_oum","cd_oumv",
                             "cid_bmv","cid_ouv","cid_oum","cid_oumv",
                             "cd_bmv_hab","cd_ouv_hab","cd_oum_hab","cd_oumv_hab",
                             "cid_bmv_hab","cid_ouv_hab","cid_oum_hab","cid_oumv_hab",
                             "cd_bmv_life","cd_ouv_life","cd_oum_life","cd_oumv_life",
                             "cid_bmv_life","cid_ouv_life","cid_oum_life","cid_oumv_life"))

models$cd_bmv_hab[2,] <- c(1,1,2,3,3)
models$cd_bmv_hab[3,] <- c(4,4,4,4,4)
models$cid_bmv_hab[2,] <- c(1,1,2,3,3,4,4,5,6,6)
models$cid_bmv_hab[3,] <- c(7,7,7,7,7,7,7,7,7,7)

models$cd_ouv_hab[2,] <- c(1,1,2,3,3)
models$cd_ouv_hab[3,] <- c(4,4,5,6,6)
models$cid_ouv_hab[2,] <- c(1,1,2,3,3,4,4,5,6,6)
models$cid_ouv_hab[3,] <- c(7,7,8,9,9,10,10,11,12,12)

models$cd_oum_hab[3,] <- c(3,3,4,5,5)
models$cid_oum_hab[3,] <- c(3,3,4,5,5,6,6,7,8,8)

models$cd_oumv_hab[2,] <- c(2,2,3,4,4)
models$cd_oumv_hab[3,] <- c(5,5,6,7,7)
models$cid_oumv_hab[2,] <- c(2,2,3,4,4,5,5,6,7,7)
models$cid_oumv_hab[3,] <- c(8,8,9,10,10,11,11,12,13,13)

models$cd_bmv_life[2,] <- c(2,1,2,2,3)
models$cd_bmv_life[3,] <- c(4,4,4,4,4)
models$cid_bmv_life[2,] <- c(2,1,2,2,3,5,4,5,5,6)
models$cid_bmv_life[3,] <- c(7,7,7,7,7,7,7,7,7,7)

models$cd_ouv_life[2,] <- c(2,1,2,2,3)
models$cd_ouv_life[3,] <- c(4,5,4,4,6)
models$cid_ouv_life[2,] <- c(2,1,2,2,3,5,4,5,5,6)
models$cid_ouv_life[3,] <- c(7,8,7,7,9,10,11,10,10,12)

models$cd_oum_life[3,] <- c(4,3,4,4,5)
models$cid_oum_life[3,] <- c(4,3,4,4,5,7,6,7,7,8)

models$cd_oumv_life[2,] <- c(3,2,3,3,4)
models$cd_oumv_life[3,] <- c(6,5,6,6,7)
models$cid_oumv_life[2,] <- c(3,2,3,3,4,6,5,6,6,7)
models$cid_oumv_life[3,] <- c(9,8,9,9,10,12,11,12,12,13)
}

# run the models (this will take a long time)
# humerus
# with sirens
c3func(HumerusSub[,c("Species","Habitat","Life","Hum_SMA")], hum.tree, models, nSim = 100, path = "hOUwie_Hum_SMA_100_StewartWiens")
c3func(HumerusSub[,c("Species","Habitat","Life","Hum_COM")], hum.tree, models, nSim = 100, path = "hOUwie_Hum_COM_100_StewartWiens")

# without sirens
c3func(HumerusSub[fem.tree$tip.label,c("Species","Habitat","Life","Hum_SMA")], fem.tree, models, nSim = 100, path = "hOUwie_Hum_SMA_100_nosiren_StewartWiens")
c3func(HumerusSub[fem.tree$tip.label,c("Species","Habitat","Life","Hum_COM")], fem.tree, models, nSim = 100, path = "hOUwie_Hum_COM_100_nosiren_StewartWiens")

# femur
c3func(FemurSub[,c("Species","Habitat","Life","Fem_SMA")], fem.tree, models, nSim = 100, path = "hOUwie_Fem_SMA_100_StewartWiens")
c3func(FemurSub[,c("Species","Habitat","Life","Fem_COM")], fem.tree, models, nSim = 100, path = "hOUwie_Fem_COM_100_StewartWiens")

#legend
# "Aq_bi"  "Aq_pd" "SAq_bi"   "T_bi"   "T_dd" 

# read in one set of models
folder = "hOUwie_Hum_SMA_100_StewartWiens"
model_set <- list()
for (i in list.files(folder, full.names = T)) {
  model_set <- c(model_set,list(readRDS(i)))
}
names(model_set) <- gsub(".RDS","",list.files(folder, full.names = F))
model_set <- model_set[names(models)]

# print model table
print(getModelTable(model_set, type = "AICc"))

# get best model(s)
best_model <- which(getModelTable(model_set, type = "AICc")[,"dAICc"] == min(getModelTable(model_set, type = "AICc")[,"dAICc"]))
best_model <- model_set[[best_model]]
best_model$solution.cont

# get the model averaged parameters
model_avg_pars <- getModelAvgParams(model_set, type = "AICc")
rownames(model_avg_pars) <- gsub(" ","_",rownames(model_avg_pars))
model_avg_pars$tip_state <- FemurSub[rownames(model_avg_pars),"HabitatLife"]
# sigma
tapply(model_avg_pars$sigma.sq,model_avg_pars$tip_state, mean)
# theta
tapply(model_avg_pars$theta,model_avg_pars$tip_state, mean)
