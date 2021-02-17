library(dplyr)
library(tidyverse)
library(stringr)
library(padr)
library(pracma)
library(cluster)
library(factoextra)


skuC5 = data.frame(read.csv("CANON_C5_EMEA.csv"))
skuPQ = data.frame(read.csv("CANON_PQ_EMEA.csv"))
skuGS = data.frame(read.csv("CANON_GS_EMEA.csv"))
sku4L = data.frame(read.csv("CANON_4L_EMEA.csv"))
sku8A = data.frame(read.csv("CANON_8A_EMEA.csv"))
sku2Q = data.frame(read.csv("CANON_2Q_EMEA.csv"))
sku2B = data.frame(read.csv("CANON_2B_EMEA.csv"))
sku1D = data.frame(read.csv("CANON_1D_EMEA.csv"))

sku2 = rbind(skuC5,skuPQ,skuGS,sku4L,sku8A,sku2Q,sku2B,sku1D)
sku3 = as.data.frame(apply(sku2[,2:ncol(sku2)], 2, function(x) {
  y = str_replace_all(x[], ",", "") #remove commas
  return(as.numeric(y)) #then convert
}))
sku3_full = cbind(sku2[,1], sku3)
names(sku3_full)[1] = "Primary.Base.Product"

sku4b1_full = sku3_full %>% filter(is.na(X19.Jul) &
                                     is.na(X19.Jun) & 
                                     is.na(X19.May) &
                                     is.na(X19.Apr) & 
                                     is.na(X19.Mar) &
                                     is.na(X19.Feb))
write.table(sku4b1_full, file="EMEA_EOL.csv",sep=",")

sku4b2_full = anti_join(sku3_full, sku4b1_full, by="Primary.Base.Product")
write.table(sku4b2_full, file="EMEA_NONEOL.csv",sep=",")


#####################################   Without PBP Info ##########################

sku4b2 = data.frame(sku4b2_full[,2:ncol(sku4b2_full)])
sku4b2 = as.data.frame(apply(sku4b2, 2, function(x) {
  y = str_replace_all(x[], ",", "") #remove commas
  return(as.numeric(y)) #then convert
}))

class(sku4b2[1,1])

sku4b2 = as.matrix(sku4b2 %>% fill_by_value(value=0))
sku4b2t =t(sku4b2)
####write.table(sku4b2t,file = "sku4trans.csv",sep=",")
####length(sku4b2t[,1])
####mean(sku4b2t[which.max(sku4b2t[,1] != 0) : length(sku4b2t[,1]),5])
####dim(sku4b2t)


sku4b3 = data.frame(sku4b2) %>% mutate(count=rowSums(.!=0))
sku4b3_full = cbind(data.frame(sku4b2_full[,1]),data.frame(sku4b3))

sku4b4_full = sku4b3_full %>% filter(count>=6)
sku4b4 = sku4b4_full[2:(ncol(sku4b4_full)-1)]
sku4b4t = t(as.matrix(sku4b4))
##################################### MEAN VAR Part ##########################

M1 = matrix(nrow = 2, ncol = ncol(sku4b4t))	
for (i in 1:ncol(sku4b4t))
{
  M1[1,i] = mean(sku4b4t[which.max(sku4b4t[,i] != 0) : length(sku4b4t[,i]),i])
  M1[2,i] = sd(sku4b4t[which.max(sku4b4t[,i] != 0) : length(sku4b4t[,i]),i])
}

class(sku4b4[1,1])
dim(M1)
rownames(M1) = c("Mean","SD")
M2 = t(M1)
####M3 = rbind(sku4b2t,M1)
M2d = data.frame(M2)
M2d = M2d %>% mutate(Cov=SD/Mean)
M2d2 = M2d %>% mutate(Covsq=Cov*Cov)
#####################################   ADI PART ##########################

M3 = matrix(nrow = 3, ncol = ncol(sku4b4t))	
for (i in 1:ncol(sku4b4t))
{
  M3[1,i] = sum(data.frame(sku4b4t[which.max(sku4b4t[,i] != 0) : length(sku4b4t[,i]),i]) !=0) 
  M3[2,i] = length((sku4b2t[which.max(sku4b4t[,i] != 0) : length(sku4b4t[,i]),i])) 
  M3[3,i] = approx_entropy(ts(sku4b4t[which.max(sku4b4t[,i] != 0) : length(sku4b4t[,i]),i],frequency=1), edim = 2, r = 0.2*sd(ts(sku4b2t[which.max(sku4b4t[,i] != 0):length(sku4b4t[,i]),i],frequency=1)), elag = 1)
}

rownames(M3) = c("NonZeroDemand","TotalDemand","Entropy")
M4= t(M3)
M4d = data.frame(M4)
M4d = M4d %>% mutate(ADI=TotalDemand/NonZeroDemand)
######################################## Pre-FINAL FILE  ###################################
M5d = cbind(M2d2,M4d)
#####################Add PBP information
M5d_full = cbind(sku4b4_full[,1],M5d)
names(M5d_full)[1] <- "Primary.Base.Product"

M6d = data.frame(cbind(M5d$Mean,M5d$Cov, M5d$ADI, M5d$Entropy))
colnames(M6d) = c("MeanVol","CoV","ADI","Entropy")
class(M6d[1,1])

######################################## CLASSIFICATION FILE ###################################
#kmeans#
fviz_nbclust(M6d,kmeans,method="wss") + geom_vline(xintercept=4,linetype=2) + ggtitle("Elbow Plot")
km.res = kmeans(M6d,4,nstart=25)
fviz_cluster(km.res,data=M6d,ggtheme=theme_minimal(),main="Cluster plot CANON EMEA")
aggregate(M6d,by=list(cluster=km.res$cluster),mean)
M6d_final = cbind(sku4b4_full[,1],M6d,km.res$cluster)
View(M6d_final)
names(M6d_final)[1] = "Primary.Base.Product"

write.table(M6d_final,file= "Canon_EMEA.csv",sep=",",row.names = FALSE)