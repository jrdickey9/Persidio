#' pop.Musse
#'
#' @description Multi-State Speciation and Extinction (MuSSE) model using bayesian based inference on parameters. Speciation and extinction parameters are informed through a population ecology and population genetics framework.
#'
#' @usage pop.Musse(phy, states, k, q.div)
#'
#' @param phy a phylogenetic tree, in ape "phylo" format.
#' @param states a numeric vector from Forest.knolls() depicting the character states of interest.
#' @param k the numeric value of character states. Three is the max number of traits for this model.
#' @param parm the parameters chosen to vary in the less constrained model. Lambda denotes speciation rates, mu for extinction rates, and qs for transition rates between stepwise states.
#'
#' @author Jonathan R. Dickey
#'
#' @references
#' Beaulieu, J. M., & O’meara, B. C. (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic biology, 65(4), 583-601.
#'
#' FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6), 1084-1092.
#'
#' Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
#'
#' Magallon S. and Sanderson M.J. 2001. Absolute diversification rates in angiospem clades. Evol. 55:1762-1780.
#'
#' Neal R.M. 2003. Slice sampling. Annals of Statistics 31:705-767.
#'
#' O'Meara, B. C., Smith, S. D., Armbruster, W. S., Harder, L. D., Hardy, C. R., Hileman, L. C., ... & Stevens, P. F. (2016). Non-equilibrium dynamics and floral trait interactions shape extant angiosperm diversity. Proc. R. Soc. B, 283(1830), 20152304.
#'
#' @import diversitree
#'
#' @importFrom graphics abline
#' @importFrom stats anova coef
#'
#' @examples
#' params<- c(.1, .15, .2, .03, .045, .06, .05, 0, .05, .05, 0, .05)
#' set.seed(5)
#' phy <- tree.musse(params, 204, x0=1)
#' K<-1000
#' p<-1/(2*K)
#' ntaxa<-204
#' gens<-100
#' r<-0.5
#' states<-tips.Dis(p,K,r,ntaxa,gens,dom="NONE")
#' k<-3
#' pop.Musse(phy,states,k,parm="lambda")
#'
#' @export

pop.Musse<-function(phy,states,k,parm=c("lambda","mu","qs")){
  if (parm=="lambda"){
    names(states)<-phy$tip.label
    phy$tip.state<-states #output of tips.Dis
    ###likelihood function and MuSSE###
    lik<-make.musse(phy,phy$tip.state,3,strict=TRUE,control=list())
    the.start<-starting.point.musse(phy,k,yule=FALSE)
    #fully constrained 'minimal' #df=3
    likb<-constrain(lik,lambda2~lambda1,lambda3~lambda1,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12) #contraining full model and allowing for stepwise transition between states
    fit.cons<-find.mle(likb,the.start[argnames(likb)])
    #now letting the speciation rates to vary #df=5
    lik.lam<-constrain(lik,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12)
    fit.lam<-find.mle(lik.lam,the.start[argnames(lik.lam)])
    #anova test
    anov.table<-anova(fit.cons,lambda=fit.lam)
    #begin mcmc
    prior<-make.prior.uniform(the.start[7],the.start[1],log=FALSE)
    #prior<-make.prior.exponential(1/(2*(the.start[1]-the.start[4]))) #make uniform prior, maybe      from object states
    samps<-mcmc(lik.lam,coef(fit.lam),nstep=1000,w=1,prior=prior,print.every=50)
    col1 <- c("red", "orange", "blue")
    profiles.plot(samples[2:4], col=col1)
    abline(v=c(lambda1,lambda2,lambda3), col=col1) #speciation rates for each character state
    the.final<-data.frame(c(lambda1,lambda2,lambda3)) #lambda parameters
    list(anov.table=anov.table, samps=samps, the.final)

  }else if (parm=="mu"){
    names(states)<-phy$tip.label
    phy$tip.state<-states #output of tips.Dis
    ###likelihood function and MuSSE###
    lik<-make.musse(phy,phy$tip.state,3,strict=TRUE,control=list())
    the.start<-starting.point.musse(phy,k,yule=FALSE)
    #fully constrained 'minimal' #df=3
    likb<-constrain(lik,lambda2~lambda1,lambda3~lambda1,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12) #contraining full model and allowing for stepwise transition between states
    fit.cons<-find.mle(likb,the.start[argnames(likb)])
    #now letting the mu rates to vary #df=5
    lik.mu<-constrain(lik,lambda2~lambda1,lambda3~lambda1,q13~0,q21~q12,q23~q12,q31~0,q32~q12)
    fit.mu<-find.mle(lik.lam,the.start[argnames(lik.lam)])
    #anova test
    anov.table<-anova(fit.cons,lambda=fit.mu)
    #begin mcmc
    prior<-make.prior.uniform(the.start[7],the.start[1],log=FALSE)
    #prior<-make.prior.exponential(1/(2*(the.start[1]-the.start[4]))) #make uniform prior, maybe      from object states
    samps<-mcmc(lik.mu,coef(fit.mu),nstep=1000,w=1,prior=prior,print.every=50)
    col1 <- c("red", "orange", "blue")
    profiles.plot(samples[2:4], col=col1)
    abline(v=c(mu1,mu2,mu3), col=col1) #speciation rates for each character state
    the.final<-data.frame(c(mu1,mu2,mu3)) #lambda parameters
    list(anov.table=anov.table, samps=samps, the.final)

  }else if (parm=="qs"){
     names(states)<-phy$tip.label
    phy$tip.state<-states #output of tips.Dis
    ###likelihood function and MuSSE###
    lik<-make.musse(phy,phy$tip.state,3,strict=TRUE,control=list())
    the.start<-starting.point.musse(phy,k,yule=FALSE)
    #fully constrained 'minimal' #df=3
    likb<-constrain(lik,lambda2~lambda1,lambda3~lambda1,mu2~mu1,mu3~mu1,q13~0,q21~q12,q23~q12,q31~0,q32~q12) #contraining full model and allowing for stepwise transition between states
    fit.cons<-find.mle(likb,the.start[argnames(likb)])
    #now letting the speciation rates to vary #df=5
    lik.q<-constrain(lik,lambda2~lambda1,lambda3~lambda1,mu2~mu1,mu3~mu1)
    fit.q<-find.mle(lik.q,the.start[argnames(lik.q)])
    #anova test
    anov.table<-anova(fit.cons,lambda=fit.q)
    #begin mcmc
    prior<-make.prior.uniform(the.start[7],the.start[1],log=FALSE)
    #prior<-make.prior.exponential(1/(2*(the.start[1]-the.start[4]))) #make uniform prior, maybe      from object states
    samps<-mcmc(lik.q,coef(fit.q),nstep=1000,w=1,prior=prior,print.every=50)
    col2 <- c("red", "orange", "blue", "purple", "black", "gray")
    profiles.plot(samples[2:4], col=col2)
    abline(v=c(q12,q13,q21,q23,q31,q32), col=col2) #speciation rates for each character state
    the.final<-data.frame(c(q12,q13,q21,q23,q31,q32)) #lambda parameters
    list(anov.table=anov.table, samps=samps, the.final)
  }
}
