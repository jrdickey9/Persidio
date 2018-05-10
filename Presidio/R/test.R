library(arules)
library(diversitree)
params<- c(.1, .15, .2, .03, .045, .06, .05, 0, .05, .05, 0, .05)
set.seed(5)
phy <- tree.musse(params, 100, x0=1)
K<-1000
p<-1/(2*K)
ntaxa<-100
gens<-100
r<-0.5
states<-tips.Dis(p,K,r,ntaxa,gens,dom="NONE")
k<-3
q.div<-5
names(states)<-phy$tip.label
phy$tip.state<-states #output of tips.Dis
###likelihood function and MuSSE###
lik<-make.musse(phy,phy$tip.state,3,strict=TRUE,control=list())
the.start<-starting.point.musse(phy,k,q.div=q.div,yule=FALSE)
#fully constrained 'minimal' #df=3
likb<-constrain(lik,lambda2~lambda1,lambda3~lambda1,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12) #contraining full model and allowing for stepwise transition between states
fit.cons<-find.mle(likb,the.start[argnames(likb)])
#now letting the speciation rates to vary #df=5
lik.lam<-constrain(lik,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12)
fit.lam<-find.mle(lik.lam,the.start[argnames(lik.lam)])
#anova test
anov.table<-anova(fit.cons,lambda=fit.lam) #some improvement
#begin mcmc
prior<-make.prior.uniform(min(the.start),max(the.start),log=FALSE) #check lower and upper bounds
#prior<-make.prior.exponential(1/(2*(the.start[1]-the.start[4]))) #make uniform prior, maybe from object states
samps<-mcmc(lik.lam,coef(fit.lam),nstep=1000,w=1,prior=prior,print.every=50)
col1 <- c("red", "orange", "blue")
profiles.plot(samps[2:4], col=col1) #ploting posterior prob of three rates
abline(v=c(lambda1,lambda2,lambda3), col=col) #speciation rates for each character state
the.final<-data.frame(c(lambda1,lambda2,lambda3)) #lambda parameters
list(anov.table=anov.table, samps=samps, datdat)
