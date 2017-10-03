library(factoextra) # segmentation pacakge
library(cluster) # segmentation pacakage
library(haven) # packages for read_sas
library(plyr)
library(dtplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(data.table)
library(pbapply)
library(FNN)
# define the work directory for your project
# equivalent to libname in SAS
setwd('/export/home/jwu/Rproject/201708_DS_Kmeans/02_data')

# verify the work directory is properly defined
getwd()

# load dataset includes all variables required for segmentation
raw <- read_sas('./analytical.sas7bdat')

# view the columns to make sure the data is correctly loaded
head(raw)

# get descriptive summary of the dataset
# equivalent as PROC MEANS  in SAS
summary(raw)

# replace missing values
impt_avg_gap <- data.frame(raw$GAP_AVG)
setnames(impt_avg_gap, "IMPT_GAP_AVG")
impt_avg_gap[is.na(impt_avg_gap)] <- 366

# cap outliers
quantile(raw$VOL_R12,c(0,0.01, 0.25, 0.5, 0.75, 0.99,1))
cap_vol_r12 <- data.frame(raw$VOL_R12)
setnames(cap_vol_r12, "CAP_VOL_R12")
cap_vol_r12[ cap_vol_r12 > 5765.504] <- 5765.504

quantile(raw$AMT_R12,c(0,0.01, 0.25, 0.5, 0.75, 0.99,1))
cap_amt_r12 <- data.frame(raw$AMT_R12)
setnames(cap_amt_r12, "CAP_AMT_R12")
cap_amt_r12[ cap_amt_r12 >  6072.16] <-  6072.16

# combine the capped and filled in missing variables to the original dataset
analytical <- cbind(raw, impt_avg_gap, cap_vol_r12, cap_amt_r12)

# standardize variables
std_cap_amt_r12 <- data.frame(scale(analytical$CAP_AMT_R12))
setnames(std_cap_amt_r12, "STD_CAP_AMT_R12")
std_cap_vol_r12 <- data.frame(scale(analytical$CAP_VOL_R12))
setnames(std_cap_vol_r12, "STD_CAP_VOL_R12")
std_impt_avg_gap <- data.frame(scale(analytical$IMPT_GAP_AVG))
setnames(std_impt_avg_gap, "STD_IMPT_GAP_AVG")

analytical <- cbind(analytical, std_cap_amt_r12, std_cap_vol_r12, std_impt_avg_gap)

# K-means segmentations
vars <- analytical[,c("STD_CAP_VOL_R12","STD_CAP_AMT_R12","STD_IMPT_GAP_AVG")]

head(vars)

# Compute and plot wss for k = 2 to k = 15 to identify optimal number of clusters
set.seed(123)
k.max <- 5 # Maximal number of clusters
data <- vars
wss <- sapply(2:k.max,
              function(k){kmeans(data, k, nstart=10 )$tot.withinss})

jpeg('NbClust.jpg')
plot(2:k.max, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
# abline(v = 3, lty =2)
dev.off()


# run the final clusters
model <- kmeans(vars,centers = 5, iter.max=1000, nstart = 5)
clus <- data.frame(model$cluster)
analytical <- cbind(analytical,clus)
head(analytical)

#check cluster size distribution and variables means of each cluster
stats <- count(analytical, "model.cluster")
stats <- cbind(stats,model$centers)
stats1 <- aggregate.data.frame(analytical[,c("CAP_AMT_R12", "CAP_VOL_R12", "IMPT_GAP_AVG")],
                               by = list(model$cluster), mean, na.rm=TRUE)
stats <- cbind(stats, stats1) 

print(stats)

# get stats for Rsquare Ratio
variable.stats<-function(v,classified){
  tot<-sd(v)
  wth<-sqrt(sum(tapply(v,classified,FUN = function (x) {sum((x-mean(x))^2)})/(length(v)-unique(classified))))
  RSq<-1-(wth/tot)^2
  Ratio<-RSq/(1-RSq)
  a<-c(tot,wth,RSq,Ratio)
  a
}

vapply(X = analytical[,c("CAP_AMT_R12", "CAP_VOL_R12", "IMPT_GAP_AVG")],FUN = variable.stats, FUN.VALUE = c(Tot.STD=0,Within.STD=0,RSQ=0,RSQRatio=0),
       classified=model$cluster)


str(analytical)

#save segments centroid for scoring next time
saveRDS(model,"./seed.RDS")

#save segment results
save(analytical,file = "./analytical_results.rda")

######################################################################
# score code going forwards
######################################################################
#load the model object saving the scoring seeds
model.km <- readRDS("./seed.RDS")

new_dataset <- analytical[,c("COLLECTOR_KEY", "STD_CAP_VOL_R12", "STD_CAP_AMT_R12", "STD_IMPT_GAP_AVG")]

summary(new_dataset)
print(model.km$centers)
library(FNN)
pred.knn <- data.frame(get.knnx(model.km$centers, new_dataset[,2:4], 5)$nn.index[,1])

count(pred.knn)

setnames(pred, "cluster_score")

pred <- cbind(new_dataset,pred)

# skip code below, just for demostration purpose to validate scoring code is working
compare <- data.frame(cbind(pred$cluster_score, analytical$model.cluster))
names(compare) <- c("cluster_score", "cluster_model")

head(compare)

# Load function
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")
crosstab(compare, row.vars = "cluster_score", col.vars = "cluster_model", type='f')
