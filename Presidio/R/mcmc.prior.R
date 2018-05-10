library(arules)
mcmc.prior<- function(K,r,gens,dom=c("NONE","AA","Aa","aa"),log=TRUE){
#creating [Q] transition matrix to inform prior distribution
    ####ex.Rescue proportion of A allele [0,(1/3)) or "1"
    p1<-runif(1,0,(1/3))
    dist1<-NA
    for(i in 1:gens){
      q.ab<-runif(1,min=min(ex.Rescue(p,K,r,gens)[[1]]),max=max(ex.Rescue(p,K,r,gens)[[1]]))
      dist1[i]<-q.ab
    }
    edges1<-c(min(dist1),quantile(dist1,(1/3)),quantile(dist1,(2/3)),max(dist1)) #creating bins going      from "1" to where? "1", "2", or "3"?
    unname(edges1)
    dat1<-discretize(dist1,edges1, method="fixed", labels=c("q11","q12","q13")) #factor
    frame1<-data.frame(dist1,dat1)
    colnames(frame1)<-c("p","bins")
    q11<-subset(frame1,bins=="q11")
    q11<-q11[,1]
    q12<-subset(frame1,bins=="q12")
    q12<-q12[,1]
    q13<-subset(frame1,bins=="q13")
    q13<-q13[,1]

    ####ex.Rescue proportion of A allele [(1/3),(2/3)) or "2"
    p2<-runif(1,(1/3),(2/3))
    dist2<-NA
    for(i in 1:gens){
      q.bc<-runif(1,min=min(ex.Rescue(p,K,r,gens)[[1]]),max=max(ex.Rescue(p,K,r,gens)[[1]]))
      dist2[i]<-q.bc
    }
    edges2<-c(min(dist2),quantile(dist2,(1/3)),quantile(dist2,(2/3)),max(dist2)) #creating bins going      from "1" to where? "1", "2", or "3"?
    unname(edges2)
    dat2<-discretize(dist2,edges2, method="fixed", labels=c("q21","q22","q23")) #factor
    frame2<-data.frame(dist2,dat2)
    colnames(frame2)<-c("p","bins")
    q21<-subset(frame2,bins=="q21")
    q21<-q21[,1]
    q22<-subset(frame2,bins=="q22")
    q22<-q22[,1]
    q23<-subset(frame2,bins=="q23")
    q23<-q23[,1]


    ####ex.Rescue proportion of A allele [(2/3),1] or "3"
    p3<-runif(1,(2/3),(1))
    dist3<-NA
    for(i in 1:gens){
      q.cd<-runif(1,min=min(ex.Rescue(p,K,r,gens)[[1]]),max=max(ex.Rescue(p,K,r,gens)[[1]]))
      dist3[i]<-q.cd
    }
    edges3<-c(min(dist3),quantile(dist3,(1/3)),quantile(dist3,(2/3)),max(dist3)) #creating bins going      from "1" to where? "1", "2", or "3"?
    unname(edges3)
    dat3<-discretize(dist3,edges3, method="fixed", labels=c("q31","q32","q33")) #factor
    frame3<-data.frame(dist3,dat3)
    colnames(frame3)<-c("p","bins")
    q31<-subset(frame3,bins=="q31")
    q31<-q31[,1]
    q32<-subset(frame3,bins=="q32")
    q32<-q32[,1]
    q33<-subset(frame3,bins=="q33")
    q33<-q33[,1]

    #qmat<-matrix(c(q11,q12,q13,q21,q22,q23,q31,q32,q33)) not sure if I should leave as distribution of possible transitions given a runif, or if I should take the mean of the distribution. Simple use this line and add mean(q33[,1]) to above lines for mean.
    upper<-c(max(q11,q12,q13,q21,q22,q23,q31,q32,q33)) #is there a way I could just do this for the whole matrix?
    lower<-c(min(q11,q12,q13,q21,q22,q23,q31,q32,q33))
    if (length(lower) == 2 && missing(upper)) {
      upper <- lower[2]
      lower <- lower[1]
    }
    n <- length(lower)
    if (length(upper) != n)
      stop("'lower' and 'upper' both be the same length")
    p.in <- 1/(upper - lower)
    p.out <- 0
    if (log) {
      p.in <- log(p.in)
      p.out <- -Inf
    }
    thisguy<-fucntion(x){ #some function here that puts out .... for mcmc prior #i dont understand x in original function what does that mean, some argument x?
      ret <- rep(p.in, length.out = length(x))
      ret[x < lower | x > upper] <- p.out
      if (log)
        sum(ret)
      else prod(ret)
    }
    return(thisguy)
 }


