##########################################################################################
### LINEAR REGRESSION with LM - JAGS & INLA (integrated nested laplace approximation) ###
##########################################################################################
#Autor: Arnau Escapa

########################
### Data Generation ###
########################
trueA<-0
trueB<-5
trueSd<-10
sampleSize<-31
# create independent x-values 
x<-(-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to a + b*x + N(0,sd)
y<-trueA+trueB*x+rnorm(n=sampleSize,mean=0,sd=trueSd)
plot(x,y,pch='+',col="blue",main="Test Data")
abline(trueA,trueB,lwd=2.5,col="cyan")

##############################
### Frequentist statistics ###
##############################
lm.1<-lm(y~x)
plot(x,y,pch='+',col="blue",main="Test Data")
abline(trueA,trueB,lwd=2.5,col="cyan")
abline(lm.1,lwd=2.5,col="blue")
legend("topleft",c("Generating line","Regression line"), lwd=2.5,col=c("cyan","blue"))
summary(lm.1)
a.hat<-as.numeric(lm.1$coefficients[1])
b.hat<-as.numeric(lm.1$coefficients[2])
lm.1.anova<-anova(lm.1)
lm.1.anova
# str(lm.1.anova)
sigma2.hat<-lm.1.anova$Sum[2]/lm.1.anova$Df[2]
sigma.hat<-sqrt(sigma2.hat)
sprintf("Estimated intercept = %f",round(a.hat,3))
sprintf("Estimated slope = %f",round(b.hat,3))
sprintf("Estimated std. deviation = %f",round(sigma.hat,3))

####################
#### JAGS model ####
####################
require(R2jags)
data = list('x' = x, 'y' = y, 'n' = sampleSize)

model.text <- "model{

#  Likelihood
for(i in 1:n){
y[i]   ~ dnorm(mu[i], pow(sd, -2))###atencion !!! JAGS works with precision. Recall presicion= 1/sigma^2
mu[i] <- a + b*x[i]
}

# Prior 
a ~ dnorm(0, pow(5, -2)) 
b ~ dunif(0, 10)
sd ~ dunif(0, 30)      
}"

model.spec<-textConnection(model.text)

jags_model <- jags.model(model.spec, data = data, n.chains = 4, n.adapt = 100)

coda_sample <- coda.samples(jags_model, 
        variable.names=c("a","b","sd"), 
        n.iter=20000, progress.bar="none")

summary(coda_sample)
plot(coda_sample)
########################
      ### INLA ###
########################
require(INLA)
imod <- inla(y ~ x , family="gaussian", data=data)
imod$summary.fixed
imod$summary.hyperpar

plot(imod)

###checking the priors paramaters
inla.set.control.fixed.default()[c('mean','prec')]

#inla - setting prior parameters
imod <- inla(y ~ x , family="gaussian", control.fixed=list(prec=1/25), data=data)
#summary(imod)
imod$summary.fixed
imod$summary.hyperpar
plot(imod) 

###################
### Predictions ###
###################
data <- data.frame(y=round(y,1), x)  #dataset
head(data)
newdata <- data.frame(x=seq(min(data$x, na.rm=TRUE), max(data$x, na.rm=TRUE), len=100), y=NA)
data.pred <- rbind(data, newdata)

imod.pred <- inla(y~x, data=data.pred, control.fixed=list(prec=1/25),control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
#extract the predictor values associated with the appended data
newdata <- cbind(newdata, imod.pred$summary.linear.predictor[(nrow(data)+1):nrow(data.pred),])
head(newdata)
#install.packages("reshape")
library(reshape)
newdata <- reshape:::rename(newdata, c("0.025quant"="lower", "0.975quant"="upper"))
library(ggplot2)
ggplot(newdata, aes(y=mean, x=x)) + geom_point(data=data, aes(y=y)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) +
  geom_line() + theme_classic()
