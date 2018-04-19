#' het.Adv
#'
#' @description het.Adv, a function describing selection on a population given a single climate change event that affects temperature by a 2 degree threshold. Exploring with the geontypes AA, Aa, and aa. This function depicts heterozygous advantage while estimating population size.
#'
#' @usage het.Adv(p,K,r,gens)
#'
#' @param p an object describing the initial frequency of A
#' @param r a number giving the population growth rate
#' @param gens is the number of generations
#' @param K describes the carrying capacity
#'
#' @author Jonathan Dickey
#'
#' @references
#' Bay, R. A., Rose, N. H., Logan, C. A., & Palumbi, S. R. (2017). Genomic models predict successful coral adaptation if future ocean warming rates are reduced. Science advances, 3(11), e1701413.
#'
#' Gomulkiewicz, R., & Holt, R. D. (1995). When does evolution by natural selection prevent extinction?. Evolution, 49(1), 201-207.
#'
#' Lindsey, H. A., Gallie, J., Taylor, S., & Kerr, B. (2013). Evolutionary rescue from extinction is contingent on a lower rate of environmental change. Nature, 494(7438), 463.
#'
#' @importFrom graphics stats
#'
#' @examples
#' K<-1000
#' p<-1/(2*K)
#' het.Adv(p,K,2,100)
#'
#' @export

het.Adv<-function(p,r,gens){
  pA<-rep(NA,gens)
  Nt<-rep(NA,gens)
  NAA<-K*p^2
  NAa<-K*2*p*q
  Naa<-K*q^2
  for(i in 1:gens){
    Temp<-(0:gens)/25
    wAA<-0.95-0.25*Temp #fitness function
    wAa<-0.75+0.02*Temp
    waa<-0.45+0.25*Temp
    ######### (d) #########
    #wAa<-wAA #playing around with the fitness function here. adding dominance to each allele.
    #wAa<-waa
    wAa<-colMeans(rbind(wAA,waa))
    #######################
    sAA<-NAA*wAA #expected number of survival
    sAa<-NAa*wAa
    saa<-Naa*waa
    Ns<-sAA+sAa+saa #current population size (those who have survived)
    Nt[i]<-Ns[i]*(1+r*(K-Ns[i])/K) #total population size
    if(Nt[i]<1){ #simulating extinction.
      break
    }
    pA[i]<-(sAA[i]+sAa[i]/2)/Ns[i] #surviving
    qa<-1-pA[i] #reproducing
    NAA<-Nt[i]*pA[i]^2
    NAa<-Nt[i]*2*pA[i]*qa #now the alleles are changing based on hardy weinburg each gen *simulating evolution*
    Naa<-Nt[i]*qa^2
  }
  list(pA=pA,Nt=Nt,plot(Nt)) #please excuse the graphs simplicity.
}
