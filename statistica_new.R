#Basic analysis of amplicon data

setwd("C:/Users/Admin/Desktop/CNR/Monticchio/Fazi/stats/stats/new/")

#Bacteria
OTUList<-read.delim("otu_bact.txt")
row.names(OTUList)<-OTUList[,1]
OTUList<-OTUList[,-1]
sort<-c("LG1","LG3","LG6","LG8", "LG10","LP1","LP3","LP6","LP8","LP10")
OTUList<-OTUList[,sort]
taxaC<-read.csv("taxa_bac_silva.csv")
row.names(taxaC)<-taxaC[,1]
taxaC<-taxaC[,-1]

taxa<-subset(taxaC, taxaC$class!="Chloroplast")

OTUList<-OTUList[row.names(OTUList)%in%row.names(taxa),] 
tOTUList<-as.data.frame(t(OTUList))
samples<-colnames(OTUList)
SumSample<-colSums(OTUList)
library("GUniFrac")
rareOTU<-Rarefy(tOTUList, depth = min(SumSample))  #Rarefaction
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)
totalSeq<-rowSums(rffOTU)
trffOTU<-t(rffOTU)
raOTU <- rffOTU[ ,colSums(rffOTU)!=0]
traOTU<-as.data.frame(t(raOTU))
paOTU<-raOTU
paOTU[paOTU>0]=1
taxa<-taxa[row.names(taxa)%in%row.names(traOTU),] 


#Archaea
ArOTUList<-read.delim("OTU_arch.txt")
row.names(ArOTUList)<-ArOTUList[,1]
ArOTUList<-ArOTUList[,-1]
ArOTUList<-ArOTUList[,sort]
taxaA<-read.csv("taxa_ar_silva.csv")
row.names(taxaA)<-taxaA[,1]
taxaA<-taxaA[,-1]
tArOTUList<-as.data.frame(t(ArOTUList))
samples<-colnames(ArOTUList)
SumSampleA<-colSums(ArOTUList)

rareArOTU<-Rarefy(tArOTUList, depth = min(SumSampleA))  #Rarefaction
rffArOTU<-as.data.frame(rareArOTU$otu.tab.rff)
totalSeq<-rowSums(rffArOTU)
trffArOTU<-t(rffArOTU)
raArOTU <- rffArOTU[ ,colSums(rffArOTU)!=0]
traArOTU<-as.data.frame(t(raArOTU))
paArOTU<-raArOTU
paArOTU[paArOTU>0]=1
taxaA<-taxaA[row.names(taxaA)%in%row.names(traArOTU),] 



variables<-read.delim("variables_new.txt")
row.names(variables)<-variables[,1]
variables<-variables[,-1]
y<-as.factor(variables$Basin)
z<-as.double(variables$depth)

#Alpha diversity


#Bacteria
alpha<-rowSums(paOTU) #Number of OTUs
alphaAr<-rowSums(paArOTU)
par(mfrow=c(2,2))
plot(y,alpha, main="number of bacteria zOTUs", xlab="basin", ylab="#zOTUs")
plot(y,alphaAr, main="number of archaea OTUs", xlab="substrate", ylab="#zOTUs")

waterbs<-glm(formula = alpha~y, family = poisson)
summary(waterbs)
waterbs<-glm(formula = alphaAr~y, family = poisson)
summary(waterbs)

plot(alpha,z, xlab="#OTUs", ylab="depth (m)", ylim=rev(range(z)), col=y,pch=16, cex=2)
plot(alphaAr, z,  xlab="#OTUs", ylab="depth (m)", ylim=rev(range(z)), col=y,pch=16, cex=2)

waterbs<-glm(formula = alpha~z, family = poisson)
summary(waterbs)
waterbs<-glm(formula = alphaAr~z, family = poisson)
summary(waterbs)

library(vegan)


#Beta diversity 
library("betapart")
par(mfrow=c(1,2))
betacore<-betapart.core(paOTU)
betamulti<-beta.multi(betacore) #the total multi-site dissimilarity across the sites, and its turnover and nestedness components
betapair<-beta.pair(paOTU) #beta.pair returns three matrices containing the pairwise between-site values of each component of beta diversity

shared<-as.data.frame(betacore$shared)
write.csv(shared, "share_bacteria.csv")
plot(hclust(betapair$beta.sor, method="complete"),hang=-1, main='Beta sorensen Bacteria', sub='', xlab='', cex=0.5) #plot cluster analysis of betapair

par(mfrow=c(2,2))
betacoreA<-betapart.core(paArOTU)
betamultiA<-beta.multi(betacoreA) #the total multi-site dissimilarity across the sites, and its turnover and nestedness components
betapairA<-beta.pair(paArOTU) #beta.pair returns three matrices containing the pairwise between-site values of each component of beta diversity
plot(hclust(betapairA$beta.sor, method="complete"),hang=-1, main='Beta sorensen Archaea', sub='', xlab='', cex=0.5) #plot cluster analysis of betapair
mantel(betapair$beta.sor,betapairA$beta.sor)

sharedA<-as.data.frame(betacoreA$shared)
write.csv(sharedA, "share_archaea.csv")

library(vegan)
#ADONIS in vegan: y~a where y is a dist matrix (dissimilarity) 
#Strata: If the experimental design has nestedness, then use strata to test hypotheses. For instance, imagine we are testing the whether a plant community is influenced by nitrate amendments, and we have two replicate plots at each of two levels of nitrate (0, 10 ppm). We have replicated the experiment in three fields with (perhaps) different average productivity. In this design, we would need to specify strata = field so that randomizations occur only within each field and not across all fields
adonis1<-adonis(betapair$beta.sor~z+y, permutations=9999)
adonis1
out<-capture.output(adonis1)
write(out, "10_adonisall.txt", append=TRUE)

adonis2<-adonis(betapairA$beta.sor~z+y, permutations=9999)
adonis2

par(mfrow=c(2,1))

betabray<-vegdist(raOTU,method="bray")
plot(hclust(betabray, method="average"), hang=-1, main='A. Bray-Curtis Bacteria', sub='', xlab='', cex=1) #plot cluster analysis of betapair
betabrayA<-vegdist(raArOTU,method="bray")
plot(hclust(betabrayA, method="complete"),hang=-1, main='B. Bray-Curtis Archaea', sub='', xlab='', cex=1) #plot cluster analysis of betapair
mantel(betabrayA,betabray)
adonis3<-adonis(betabray~z+y, permutations=9999, strata=variables$municipality)
adonis4<-adonis(betabrayA~z+y, permutations=9999, strata=variables$municipality)


plot(cldw,cldwA)



library(ggplot2)

#Taxonomy all bacterial otus###################

###TAXONOMY
library("reshape2")

#taxa as expample on family level only
aggT<-do.call("rbind", as.list(by(traOTU[,], taxa$family, colSums))) #make the sum of the OTU table based on the taxa$f column meaning the family
alldataT<-as.data.frame(rowSums(aggT)) #calculated the total abundace of each family


rareST4<-subset(aggT, alldataT<2000) #make a dataframe with less abundant families
abundST4<-subset(aggT, alldataT>=2000) #...and abundant families
sumrareST4<-as.data.frame(colSums(rareST4)) # sum of rare familie per sample
sumrareST4<-t(sumrareST4)
row.names(sumrareST4)<-"other"  #name them other
colnames(sumrareST4)<-colnames(abundST4)
allST4<-as.data.frame(rbind(abundST4, sumrareST4)) 
rownames(allST4)[rownames(allST4) == ""] <- "unassigned"

#Plot as heat map
mST<-as.matrix(allST4) #trasfer to matrix
mfam<-melt(mST)
#ggplot(mfam,aes(x = Var2, y = Var1 ,size = value)) + geom_point(alpha=0.4) +scale_size() +theme_light()

c22 <- c(
  "dodgerblue2",  # red
  "green4",
  "plum4", # purple
  "orange2", # orange
  "mistyrose1",
  "snow1", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon",  "steelblue4",
  "darkturquoise", "green1", "yellow4", 
   '#fabebe',	"yellow3",'#ffd8b1',	"peachpuff3","darkorange4",'#000075'
)
taxa_bac<-ggplot(mfam,aes(x = Var2, y = value ,fill = Var1)) + 
  geom_bar(position="fill", stat="identity") + 
  guides(fill=guide_legend(ncol=1, reverse = F, title="Bacteria family")) + 
  scale_x_discrete(name="") + 
  scale_y_continuous(name="Relative Abundance") + 
  theme(axis.title.y=element_text(colour="#999999"))+
  scale_fill_manual(values = c22)


#####
#Archaea

#taxa as expample on family level only
tArOTU
aggT<-do.call("rbind", as.list(by(traArOTU[,], taxaA$family, colSums))) #make the sum of the OTU table based on the taxa$f column meaning the family
alldataT<-as.data.frame(rowSums(aggT)) #calculated the total abundace of each family


rareST4<-subset(aggT, alldataT<1000) #make a dataframe with less abundant families
abundST4<-subset(aggT, alldataT>=1000) #...and abundant families
sumrareST4<-as.data.frame(colSums(rareST4)) # sum of rare familie per sample
sumrareST4<-t(sumrareST4)
row.names(sumrareST4)<-"other"  #name them other
colnames(sumrareST4)<-colnames(abundST4)
allST4<-as.data.frame(rbind(abundST4, sumrareST4)) 
rownames(allST4)[rownames(allST4) == ""] <- "unassigned"
mST<-as.matrix(allST4) #trasfer to matrix
mfam<-melt(mST)
ggplot(mfam,aes(x = Var2, y = Var1 ,size = value)) + geom_point(alpha=0.4) +scale_size() +theme_light()

taxa_arc<-ggplot(mfam,aes(x = Var2, y = value ,fill = Var1)) + 
  geom_bar(position="fill", stat="identity") + 
  guides(fill=guide_legend(ncol=1, reverse = F, title="Archaea family")) + 
  scale_x_discrete(name="") + 
  scale_y_continuous(name="Relative Abundance") + 
  theme(axis.title.y=element_text(colour="#999999"))+
  scale_fill_manual(values = c22)


library(gridExtra)
library(cowplot)

theme_set(theme_minimal())



plot_grid(
  taxa_bac+ theme(legend.justification = c(0,1))
  , taxa_arc+ theme(legend.justification = c(0,1))
  , align = "v"
  , ncol = 1
  ,labels = c('A', 'B')
)


#genus

aggT<-do.call("rbind", as.list(by(traOTU[,], taxa$genus, colSums))) #make the sum of the OTU table based on the taxa$f column meaning the family
alldataT<-as.data.frame(rowSums(aggT)) #calculated the total abundace of each family

rareST4<-subset(aggT, alldataT<1000) #make a dataframe with less abundant families
abundST4<-subset(aggT, alldataT>=1000) #...and abundant families
sumrareST4<-as.data.frame(colSums(rareST4)) # sum of rare familie per sample
sumrareST4<-t(sumrareST4)
row.names(sumrareST4)<-"other"  #name them other
colnames(sumrareST4)<-colnames(abundST4)
allST4<-as.data.frame(rbind(abundST4, sumrareST4)) 
rownames(allST4)[rownames(allST4) == ""] <- "unassigned"
#Plot as heat map
tallST4B<-t(allST4)
mst4<-as.matrix(allST4)
mtalB<-melt(mst4)

BG<-ggplot(mtalB,aes(x = Var2, y = value,fill = Var1)) + 
  geom_bar(position="fill", stat="identity") + 
  guides(fill=guide_legend(ncol=1, reverse = F, title="Bacteria Genus")) + 
  scale_x_discrete(name="") + 
  scale_y_continuous(name="Relative Abundance") + 
  theme(axis.title.y=element_text(colour="#999999"))+
  scale_fill_manual(values = c22)


aggT<-do.call("rbind", as.list(by(traArOTU[,], taxaA$genus, colSums))) #make the sum of the OTU table based on the taxa$f column meaning the family
alldataT<-as.data.frame(rowSums(aggT)) #calculated the total abundace of each family

rareST4<-subset(aggT, alldataT<700) #make a dataframe with less abundant families
abundST4<-subset(aggT, alldataT>=700) #...and abundant families
sumrareST4<-as.data.frame(colSums(rareST4)) # sum of rare familie per sample
sumrareST4<-t(sumrareST4)
row.names(sumrareST4)<-"other"  #name them other
colnames(sumrareST4)<-colnames(abundST4)
allST4<-as.data.frame(rbind(abundST4, sumrareST4)) 
rownames(allST4)[rownames(allST4) == ""] <- "unassigned"
tallST4A<-t(allST4)
mst4<-as.matrix(allST4)
mtalA<-melt(mst4)

AG<-ggplot(mtalA,aes(x = Var2, y = value,fill = Var1)) + 
  geom_bar(position="fill", stat="identity") + 
  guides(fill=guide_legend(ncol=1, reverse = F, title="Archaea Genus")) + 
  scale_x_discrete(name="") + 
  scale_y_continuous(name="Relative Abundance") + 
  theme(axis.title.y=element_text(colour="#999999"))+
  scale_fill_manual(values = c22)

plot_grid(
  BG+ theme(legend.justification = c(0,1))
  , AG+ theme(legend.justification = c(0,1))
  , align = "v"
  , ncol = 1
  ,labels = c('A', 'B')
)

#Genera with relation to depth
corA<-cbind(tallST4A, z)
corz<-as.data.frame(cor(corA))
corArch<-as.data.frame(corz$z)
rownames(corArch)<-rownames(corz)

corB<-cbind(tallST4B, z)
corz<-as.data.frame(cor(corB))
corBac<-as.data.frame(corz$z)
rownames(corBac)<-rownames(corz)

tB<-as.data.frame(tallST4B)
tA<-as.data.frame(tallST4A)
 
library("dplyr")
par(mfrow=c(2,4))
for (i in 1:14) {
plot(tB[,i],z, ylab="depth", xlab="reads", main=colnames(tB)[i], ylim=rev(range(z)), col=y,  pch=16, cex=3)
}
dev.off()
par(mfrow=c(2,4))
for (i in 1:15) {
  plot(tA[,i],z, ylab="depth", xlab="reads", main=colnames(tA)[i], ylim=rev(range(z)), col=y,  pch=16, cex=3)
}  
plot(tB$Candidatus_Planktophila,z, ylab="depth", xlab="reads", main="Candidatus Planktophila /Actinobacteria", ylim=rev(range(z)), col=y,  pch=16)
plot(tB$`Christensenellaceae_R-7_group`,z, ylab="depth", xlab="reads", main="Christensenellaceae / Clostridia", ylim=rev(range(z)), col=y,  pch=16)
plot(tA$Methanoculleus,z, ylab="depth", xlab="reads", main="Methanoculleus", ylim=rev(range(z)), col=y,  pch=16)
plot(tA$Methanomethylovorans,z, ylab="depth", xlab="reads", main="Methanomethylovorans", ylim=rev(range(z)), col=y,  pch=16)
plot(tA$Methanothermobacter,z, ylab="depth", xlab="reads", main="Methanothermobacter", ylim=rev(range(z)), col=y,  pch=16)

###Chloroplasts

taxaClo<-subset(taxaC, taxaC$class=="Chloroplast")

OTUListC<-OTUList[row.names(OTUList)%in%row.names(taxaClo),] 
Chl<-colSums(OTUList)
plot(Chl,variables$depth, cex=2, xlab="Chloroplast reads", ylab="depth (m)", pch=16,ylim=rev(range(z)), col=y)
        