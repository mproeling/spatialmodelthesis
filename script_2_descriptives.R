rm(list = ls(all = TRUE))

library(igraph)
library(spatialprobit)
library(corrplot)

# this model follows after script_createfeatures.R

setwd("C:\\Users\\mark\\Downloads\\ToreOpsahl")

FBattributes = read.table("attributes.csv", sep = ";", h = T)
FBattributes$gender01 = FBattributes$gender
FBattributes[FBattributes$gender01 == 1, ]$gender01 <- 0 # recode gender
FBattributes[FBattributes$gender01 == 2, ]$gender01 <- 1 # recode gender

FBattributes = FBattributes[order(FBattributes$id), ] # sort on ID to match W
# add intercept, remove ID and send.friends (equal to outdegree) 
FBattributes = as.data.frame(cbind(1, FBattributes[,-c(1, which(names(FBattributes) == "send.friends"))])) 
colnames(FBattributes)[1] = "Intercept"
head(FBattributes)

# correlation matrix without intercept and gender
FB.cor.matrix = cor(FBattributes[,-c(1, which(names(FBattributes) == "gender"),
                                        which(names(FBattributes) == "gender01"))])

# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
corrplot(FB.cor.matrix, method = "number")

# test mean differences
mean(FBattributes[FBattributes$gender == 1, ]$yearofstudy)
#[1] 2.495528
mean(FBattributes[FBattributes$gender == 2, ]$yearofstudy)
#[1] 2.238156
t.test(FBattributes$yearofstudy ~ FBattributes$gender)
t.test(FBattributes$outdegree ~ FBattributes$gender)
t.test(FBattributes$indegree ~ FBattributes$gender)
t.test(FBattributes$popular ~ FBattributes$gender)
t.test(FBattributes$averageChar ~ FBattributes$gender)
t.test(FBattributes$dayactive ~ FBattributes$gender)
t.test(FBattributes$send.day34th.friends ~ FBattributes$gender)

# check 
isSymmetric(W)
fit.all <- sarprobit(gender ~ yearofstudy + popular + averageChar + dayactive + send.day34th.friends + send.friends,
                     W=W, data=FBattributes, ndraw=1000, burn.in = 100, showProgress=TRUE)

summary(fit.all)
