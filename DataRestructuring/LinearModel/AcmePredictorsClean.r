
install.packages("corrplot")
install.packages("Hmisc")
install.packages("PerformanceAnalytics")

library(ggplot2)
library(scales)
library(RColorBrewer)
options(digits=2)

acme <- read.csv(file="acmeupdate1.csv",head=TRUE,sep=",")
acme
colnames(acme)

vars <- cbind (acme$Branch, acme$Actions, acme$NewCommits, log(acme$Automations + 1), log(acme$Insertions), log(acme$Deletions), log(acme$FilesModified), log(acme$FilesDeleted + 1), log(acme$FilesAdded + 1), log(acme$AutomationDifficulty + 1), log(acme$ManualAddDifficulty+1), log(acme$AddDelEqMoved + 1), log(acme$ManualLineEstimate + 1), acme$ManualLineProp, log(acme$Dependency_UserIssue+1), acme$Difficulty);
dimnames(vars)[[2]] <- c("Branch", "Actions", "NewCommits", "Automations", "Insertions", "Deletions", "FMod", "FDel", "FAdd", "AutoDifficulty", "ManualAddDifficulty", "ExplainedByMove", "ManualLineEstimate", "MLProportion","UserIssueDependency", "Difficulty")

data <- data.frame(vars)

summary(data)


cor(data)

cor(data, method="spearman")

# Define a function
hiCor <- function(x, level){
  res <- cor(x,method="spearman");
  res1 <- res; res1[res<0] <- -res[res < 0];
  for (i in 1:dim(x)[2]){
    res1[i,i] <- 0;
  }
  sel <- apply(res1,1,max) > level;
  res[sel,sel];
}
hiCor(data,.8)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
library(Hmisc)
res2<-rcorr(vars)
flattenCorrMatrix(res2$r, res2$P)

library(corrplot)
res <- cor(data)
corrplot(res, type = "full", order = "hclust", 
         tl.col = "black", tl.srt = 45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)

#PCA all variables
plot(1:12,cumsum(prcomp(vars, retx=F,scale=T)$sdev^2)/sum(prcomp(vars, retx=F,scale=T)$sdev^2),ylim=c(0,1),xlab="Number of components",ylab="Fraction of variance");
res<-prcomp(vars, retx=F,scale=T)$rotation[,1:5];
resAbs <- res;
resAbs[res<0] <- -res[res<0];
for (i in 1:5)
  print(t(res[resAbs[,i]>.3,i,drop=FALSE]));

#PCA reduced set
data2 <- subset(data, select = c("Branch", "Automations", "Insertions", "FAdd", "AutoDifficulty", "ManualAddDifficulty", "UserIssueDependency", "Difficulty", "FDel"))

plot(1:9,cumsum(prcomp(data2, retx=F,scale=T)$sdev^2)/sum(prcomp(data2, retx=F,scale=T)$sdev^2),ylim=c(0,1),xlab="Number of components",ylab="Fraction of variance");
res<-prcomp(data2, retx=F,scale=T)$rotation[,1:5];
resAbs <- res;
resAbs[res<0] <- -res[res<0];
for (i in 1:5)
  print(t(res[resAbs[,i]>.3,i,drop=FALSE]));


res <- c();
vnam <- names(data2);
for (i in 2:dim(data2)[2]){
  fmla <- as.formula(paste(vnam[i],paste(vnam[-c(1,i)],collapse="+"),sep="~"));
  res <- rbind(res,c(i,round(summary(lm(fmla,data=data))$r.squared,2)));
}
row.names(res) <- vnam[res[,1]];
res[order(-res[,2]),];

################LINEAR MODEL
fmla ~ Difficulty ~ Automations+Insertions+AutoDifficulty+UserIssueDependency+FAdd+ManualAddDifficulty

#LM
mod4 <- lm(Difficulty ~ Insertions+AutoDifficulty+UserIssueDependency + FAdd,data=data1);
summary(mod4); 

anova(mod4, test="Chi");

mod5 <- lm(Difficulty ~ Automations+FAdd+AutoDifficulty+UserIssueDependency,data=data);
summary(mod5); 
drop1(mod5)

anova(mod5, test="Chi");

#Alternative sums of squares for ANOVA
#Based on residuals of remaining predictors
drop1(mod4, test="Chi");

#Variance Inflation Factor
library(car)
vif(mod4);





#################Scratch Work Below

#cols <- c("Actions","NewCommits")
#acme[cols] <- scale(acme[cols])
#acme[cols]
normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
}

acmeNorm <- as.data.frame(lapply(acme, normalize))

acmeNorm
corAN = cor(acmeNorm)
corAN

#PCA based on correlation matrix
# Pricipal Components Analysis
# entering raw data and extracting PCs 
# from the correlation matrix 
fit <- princomp(corAN, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

colnames(acmeNorm)
mod <- glm(Difficulty ~ Automations+Insertions+AutomationDifficulty+Dependency_UserIssue+FilesAdded+ManualAddDifficulty,family=gaussian,data=acmeNorm);
summary(mod); 

mod <- glm(Difficulty ~ Automations+Insertions+AutomationDifficulty+Dependency_UserIssue,family=gaussian,data=acmeNorm);
summary(mod); 





acmeNorm.pca <- prcomp(acmeNorm,
                 center = TRUE,
                 scale. = TRUE) 
summary(acmeNorm.pca)
#1st row, standard deviation associated with each PC
#2nd row, proportion of variance in data explained by each component 
#3rd row, cumulative proportion of explained variance

#PC1,PC2, PC3, PC4

predict(acmeNorm.pca, 
        newdata=tail(acmeNorm, 2))

myPCA <- prcomp(acmeNorm, scale. = F, center = F)
myPCA

sub1 <- subset(acmeNorm, select = c("Branch", "NewCommits", "Automations", "Insertions", "Deletions", "FilesModified", "AutomationDifficulty", "ManualAddDifficulty", "ManualLineEstimate", "Dependency_UserIssue", "Difficulty"))
sub1
#cor(sub1)

sub2 <- subset(acmeNorm, select = c("Branch", "Insertions", "Deletions", "FilesModified", "AutomationDifficulty", "ManualAddDifficulty", "ManualLineEstimate", "Dependency_UserIssue", "Difficulty"))
#sub2
subset.mod <- lm(Difficulty ~ Branch + Insertions + Deletions, data=sub2)
summary(subset.mod)
plot(subset.mod, which = c(1,2))

#Formula call -> predictor and target response variables with the data being used
#

########### PCA with subset 
corSub2 = cor(sub2)
fit <- princomp(corSub2, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

sub2.pca <- prcomp(sub2,
                 center = TRUE,
                 scale. = TRUE) 
summary(sub2.pca)

myPCA2 <- prcomp(acmeNorm, scale. = F, center = F)
summary(myPCA2)







#Add variables
#use of variable before hand, indicator has it been used as
#missign instances where it should have been changed

# Summary
# Log transform
# identify skewed data
# principle components
# look at correlations
# subset of best fit to model
# max of 3 or 4
# look at model
# colnames(acme)
library(PerformanceAnalytics)
options(warn=-1)
chart.Correlation(acmeNorm[,2:9],col=acmeNorm$Difficulty)

chart.Correlation(acmeNorm[,10:16],col=acmeNorm$Difficulty)

#max(acme$Actions)
#min(acme$Actions)
#length(acme$Actions)
#plot(density(acme$Actions))
#new2 <- normalized(acme)
#new1 <- scale(new2)
#cols <- c("Actions","NewCommits")
#acme[cols] <- scale(acme[cols])
#acme[cols]

sub <- subset(acme, select = c("Branch", "Actions", "NewCommits", "Automations", "Insertions", "Deletions", "FilesModified", "FilesDeleted", "FilesAdded", "AddDifficulty", "ManualLineEstimate", "Difficulty"))
sub
cor(sub)

###############OLD SCRATCH WORK#################################################




subset <- subset(acme, select = c("Branch","Automations", "Insertions", "Deletions", "FilesModified", "AddDifficulty", "ManualLineEstimate", "Difficulty"))
subset
subset.mod <- lm(Difficulty ~ AddDifficulty + Insertions + Automations, data=acme)
summary(subset.mod)
plot(subset.mod, which = c(1,2))

subset.mod2 <- lm(Difficulty ~ Automations + Insertions + Deletions + FilesModified + AddDifficulty + ManualLineEstimate, data=acme)
summary(subset.mod2)
plot(subset.mod2, which = c(1,2))

sub1 <- subset(acme, select = c("Branch", "Actions", "NewCommits", "Automations", "Difficulty"))
sub1
cor(sub1)

sub2 <- subset(acme, select = c("Branch", "Insertions", "Deletions", "FilesModified", "FilesDeleted", "FilesAdded", "Difficulty"))
sub2
cor(sub2)

plot(sub2)

sub2.mod <- lm(Difficulty ~ Insertions, data=acme)
summary(sub2.mod)

plot(sub2.mod, which = c(1,2))

sub2.mod2 <- lm(Difficulty ~ Insertions + Deletions + FilesModified, data=acme)
summary(sub2.mod2)
plot(sub2.mod2, which = c(1,2))

matplot(acme[, 1], acme[, -1])


model<-lm(acme$Difficulty~acme$Actions + acme$NewCommits + acme$Automations)
plot(model)

create_model <- function(acme,target) {
set.seed(120)
myglm <- glm(target ~ . , data=trainData, family = "binomial")
return(myglm) }

regmodel <- lm(Difficulty ~ Actions, data = acme)
regmodel
plot(regmodel)

plot(acme$Difficulty,acme$ManualLineProportion)
cor(acme$Difficulty,acme$ManualLineProportion)




