euler<-function(T=100, N=10, start, delta=0.1 ,sigma=0.1, mu=0.1, pb=FALSE){
  ##constructing times
  delta.time<-T/N
  times<-matrix(data = seq(from = 0,by = delta.time,length.out = N), nrow = 1, ncol = N)
  ##Initialisation of compartments
  nbCompartments<-length(start)
  ##initialization of the output
  out<-matrix(nrow = length(times), ncol = nbCompartments+1)
  out[1,]<-c(0,start)
  ##sigma to matrix
  if(is.matrix(sigma)) sigma.mat<-sigma else sigma.mat<-matrix(data = sigma,nrow = nbCompartments, ncol = nbCompartments)
  ##Progress bar
  if(pb==TRUE) pbPrint <- txtProgressBar(min = 0, max = length(times[1,]), style = 3)
  ##loop over time
  for(i in seq_len(length(times[1,])-1)){
    # update progress bar
    if(pb==TRUE){setTxtProgressBar(pbPrint, i)}
    ##Generating noise increments
    gam<-apply(X = sigma.mat,MARGIN = c(1,2),FUN = function(x) rgamma(n = 1,scale =x^2 ,shape = delta/x^2))
    ##generate probabilities & increments
    #initialization
    proba<-matrix(nrow = nbCompartments, ncol = nbCompartments)
    var.N<-matrix(nrow = nbCompartments, ncol = nbCompartments)
    ##Probabilities
    proba<-(1-exp(-rowSums(mu * gam)))*(mu * gam)/(rowSums(mu * gam))
    ##Normalization
    diag(proba) <- 1 - rowSums(proba)
    ##generate process increments
    ###hoehle: I guess this can't be vectorized s.t. rmultinom(n=nCompartments, size=X[i,], proba) (I don't think so buch check)
    for(c in 1:nbCompartments){
      var.N[c,]<-rmultinom(n=1, size = out[i,c+1], prob = proba[c,])
    }
    out[i+1,-1]<-colSums(var.N)
  }
  #Output: time and compartments number at each time (# of compartments +1)
  out[,1] <- times
  return(out)
}

ret <- euler(start=c(1000,20),N=1000,T=200,delta = .1)
plot(ret[,3])

