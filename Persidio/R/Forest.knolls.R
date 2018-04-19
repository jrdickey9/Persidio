#' Forest.knolls
#'
#' @description Forest.knolls, a function that estimates the proportion of the A allele in a population given a single generation of interest and n number of taxa while discreteizing the proportion of the A allele into three traits 1 (high), 2 (medium), 3 (low) which will be an indicator for population size.
#'
#' @usage Forest.knolls(p,K,r,ntaxa, gens)
#'
#' @param p a numeric object with the intended proportion of the A allele in the populationche
#' @param K a numeric object describing the carrying capacity.
#' @param r a numeric object establishing the growth rate of a population.
#' @param ntaxa a numeric object with the number of taxa at the tips.
#' @param gens a numeric object with the number of generations of interest or time in gens at the tip.
#'
#' @author Jonathan R. Dickey
#'
#' @references
#' FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6), 1084-1092.
#'
#' Ng, J., & Smith, S. D. (2014). How traits shape trees: new approaches for detecting character state‐dependent lineage diversification. Journal of evolutionary biology, 27(10), 2035-2045.
#'
#' O'Meara, B. C., Smith, S. D., Armbruster, W. S., Harder, L. D., Hardy, C. R., Hileman, L. C., ... & Stevens, P. F. (2016). Non-equilibrium dynamics and floral trait interactions shape extant angiosperm diversity. Proc. R. Soc. B, 283(1830), 20152304.
#'
#' @importFrom diversitree, ape, arules, graphics, stats
#'
#' @examples
#' K<-1000
#' p<-1/(2*K)
#' Forest.knolls(p,K,0.5,204,100)
#'
#' @export

Forest.knolls<-function(p,K,r,ntaxa,gens){
  pAmean<-mean(Twin.peaks(p,K,r,gens)$pA)
  pAsims<-replicate(ntaxa,Twin.peaks(runif(1,pAmean,min=0),K,r,gens)$pA[gens]) #maybe I should make this match
  edges<-c(min(pAsims),quantile(pAsims,(1/3)),quantile(pAsims,(2/3)),max(pAsims)) #begin discretizing, creating bins
  unname(edges)
  dat<-discretize(pAsims,edges, method="fixed", labels=c("1","2","3")) #factor
  dat<-as.numeric(levels(dat)[dat]) #numeric
  return(dat)
}
