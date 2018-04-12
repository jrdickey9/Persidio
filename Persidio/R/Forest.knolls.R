#' Forest.knolls
#' @description Forest.knolls, a function that estimates tip population size given a random uniform distribution of A allele at a single time for n number of taxa.
#' @usage Forest.knolls(p,K,r,ntaxa)
#' @param K describes the carrying capacity
#' @param r established growth rate
#' @param ntaxa the number of taxa
#' @param gens the number of generations of interest or time in gens at the tip.
#'
#' @examples
#'Forest.knolls(K,0.5,204,100)
Forest.knolls<-function(p,K,r,ntaxa,gens){
  datmean<-mean(Twin.peaks(p,K,r,gens)$pA)
  sims<-replicate(ntaxa,Twin.peaks(runif(1,datmean,min=0),K,r,gens)$Nt[gens])
}

