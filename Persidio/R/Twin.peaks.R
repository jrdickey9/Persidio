#' Twin.peaks
#' @description Twin.peaks, a function describing selection, survival and reproduction, on a population given a single climate change event that affects temperature in an island ecosystem, 2 degree threshold. Exploring with the geontypes AA, Aa, and aa. No allele has dominance and the population is in hardy weinberg equilibrium to denote genotype frequency.
#' @usage Twin.peaks(p,K,r,gens)
#' @param p an object describing the initial frequency of A
#' @param r a number giving the population growth rate
#' @param gens is the number of generations
#' @param K describes the carrying capacity
#'
#' @examples
#' Twin.peaks(p,K,2,100)
K<-1000
p<-1/(2*K)
Twin.peaks<-function(p,K,r,gens){
  K<-K
  q<-1-p
  pA<-rep(NA,gens) #numeric(length=t)
  Nt<-rep(NA,gens) #numeric (length=t)
  NAA<-K*p^2
  NAa<-K*2*p*q
  Naa<-K*q^2
  for(i in 1:gens){
    Temp<-(0:gens)/25
    wAA<-0.95-0.25*Temp #fitness function
    wAa<-0.75+0.02*Temp
    waa<-0.45+0.25*Temp

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
    NAA<-Nt[i]*pA[i]^2 #HW equilibrium
    NAa<-Nt[i]*2*pA[i]*qa #now the alleles are changing based on hardy weinburg each gen *simulating evolution*
    Naa<-Nt[i]*qa^2
  }
  list(pA=pA,Nt=Nt,plot(Nt)) #when I try adding fancy arguments in the plot function within the output some weird errors occur, please excuse the graphs simplicity.
}
