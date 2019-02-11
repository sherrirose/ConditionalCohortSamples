################################################################
####   Demonstration of observation weighting technique     ####
####   from the paper "Consistent Estimation of Propensity  ####
####   Score Functions with Oversampled Exposed Subjects"   ####
####   by Sherri Rose arxiv.org/abs/1805.07684              ####
################################################################


library(SuperLearner)
SL.randomForest.250<- function(..., size=250){SL.randomForest(...,size=size)}
SL.nnet.3<- function(..., size=3){SL.nnet(...,size=size)}
SL.nnet.5<- function(..., size=5){SL.nnet(...,size=size)}
SL.lib <- c("SL.glm","SL.randomForest.250","SL.rpart",
            "SL.nnet", "SL.nnet.3","SL.nnet.5", "SL.glmnet")

####Creating Source Population####
set.seed(980)
Npop<-50000 #number of individuals in the population
pop<-data.frame(X1 = rbinom(Npop, 1, .6), X2 = rbinom(Npop, 1, .4),
X3 = rbinom(Npop, 1, .4), X4 = rbinom(Npop, 1, .5), X5 = rbinom(Npop, 1, .4),
X6 = rbinom(Npop, 1, .5))
tru=1/(1+exp(-(3*pop$X1+1.1*pop$X2+2.2*pop$X3-1.7*pop$X4-4.8*pop$X5-3.7*pop$X6)))
pop<-transform(pop, E = rbinom(Npop, 1, tru), tru=tru)
npop<-sum(pop$E)  #number of exposed in the population
nCpop<-Npop-npop  #number of unexposed controls in the population
w<-npop/Npop #probability of exposure in the population
unexposed.pop<-data.frame(matrix(data=NA, nrow=nCpop, ncol=dim(pop)[2]))
exposed.pop<-data.frame(matrix(data=NA, nrow=npop, ncol=dim(pop)[2]))
items<-c("X1", "X2", "X3", "X4", "X5", "X6","E", "tru")
colnames(unexposed.pop)<-items;colnames(exposed.pop)<-items

####Sorting the Population into Exposed and Unexposed Populations####
pop.sort<-pop[order(pop$E),]
for(i in 1:Npop){
if (pop.sort[i,"E"]==0) {unexposed.pop[i,]=pop.sort[i,]}
}
pop.sort<-pop[order(pop$E, decreasing=TRUE),]
for(i in 1:Npop){
if (pop.sort[i,"E"]==1) {exposed.pop[i,]=pop.sort[i,]}
}

####Setting Values to Create Samples####
n<-200 #number of exposed subjects
C<-1 #number of controls matched to each exposed subject
nC<-n*C #number of controls in sample
N<-n*(C+1) #number of observations in study sample

####Creating Sample Conditional on Exposure####
exposed.sample<-exposed.pop[(sample((1:npop), n, replace=FALSE)),]
unexposed.sample<-unexposed.pop[(sample((1:nCpop),nC, replace=FALSE)),]
conditionalsample<-rbind(exposed.sample, unexposed.sample)

####Creating Random Sample####
randomsample<-pop[(sample((1:Npop), N, replace=FALSE)),]

####Observation Weights####
wC<-(1-w)/C
weight<-ifelse(conditionalsample$E==1,w, wC)

####Ensemble Estimation####
SL.random <- SuperLearner(Y=randomsample[,7], X=randomsample[,1:6],
                           family="binomial", SL.library=SL.lib,
                           method="method.NNLS", verbose=T)
SL.weights <- SuperLearner(Y=conditionalsample[,7], X=conditionalsample[,1:6],
                           family="binomial", SL.library=SL.lib,
                           method="method.NNLS", verbose=T, obsWeights=weight)

SL.noweights <- SuperLearner(Y=conditionalsample[,7], X=conditionalsample[,1:6],
                            family="binomial", SL.library=SL.lib,
                            method="method.NNLS", verbose=T)

####Results: Bias###
mean(randomsample$tru-SL.random$SL.predict)
mean(conditionalsample$tru-SL.weights$SL.predict)
mean(conditionalsample$tru-SL.noweights$SL.predict)

####Results: MSE###
mean((randomsample$tru-SL.random$SL.predict)^2)
mean((conditionalsample$tru-SL.weights$SL.predict)^2)
mean((conditionalsample$tru-SL.noweights$SL.predict)^2)

####Results: Algorithm Coefficients###
SL.random$coef
SL.weights$coef
SL.noweights$coef
