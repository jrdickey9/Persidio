#' Dolores
#' @descrtiption Biogeography model using maximum likelihood based inference on parameters and ancestral states. Speciation and extinction parameters are informed through a population ecology framework. In this biogeographic context, a climate event has occured that shifts temperature 2 degrees.
#' @usage Dolores(phy, dat, tp, )
#' @param phy a phylogenetic tree, in ape "phylo" format
#' @param dat a data matrix contatining species information
#' @param tp a list produced by Twin.peaks() with pA [1] and Nt [2]
#' @author Jonathan R. Dickey
#' @references
#' Beaulieu, J. M., & O’meara, B. C. (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic biology, 65(4), 583-601.
#'
#' Goldberg E.E., Lancaster L.T., and Ree R.H. 2011. Phylogenetic inference of reciprocal effects between geographic range evolution and diversification. Syst. Biol. 60:451-465.
#'
#' Wayne P. Maddison, Peter E. Midford, Sarah P. Otto; Estimating a Binary Character's Effect on Speciation and Extinction, Systematic Biology, Volume 56, Issue 5, 1 October 2007, Pages 701–710.
#'
#' Magallon S. and Sanderson M.J. 2001. Absolute diversification rates in angiospem clades. Evol. 55:1762-1780.
#' @example Dolores()
Dolores<-function(phy, dat, tp, ){

}
