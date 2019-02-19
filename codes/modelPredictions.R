volume = function(t, ag , am, bg, bm, th, pdef=0) { (ag/th * (1 - (th*exp(-bg*t) - bg*exp(-th*t))/(th - bg)) - am/th*(1 - (th*exp(-bm*t) - bm*exp(-th*t))/(th - bm)))*(1-pdef) }

deltaVPrediction <- function(V0, Vext, om0, rho,e){
  # get ome (proportion of commercial volume in volume loss due to logging operations)
  mu = om0^(1-rho); Var = (1-mu)*mu^2*e; 
  alpha =   (mu/Var - 1/mu) * mu^2;beta = alpha * (1/mu - 1);
  ome = apply(cbind(alpha,beta,om0), 1, function(X) {
    if (X[3]>0 & X[3]<1) rtrunc(1,"beta",shape1=X[1],shape2=X[2],a=X[3]) else if (X[3]>1) 1 else 0
  })
  
  deltaV = apply(cbind(Vext/ome,V0), 1, min)
  deltaV[om0==0] <- 0
  return(deltaV)
}

t0Prediction <- function(ti, deltaV, aG,aM,bP,bM,theta,pdef=0){
  equ = function(x) volume(x,aG,aM,bP,bM,theta,pdef) - (volume(ti,aG,aM,bP,bM,theta,pdef) - deltaV)
  return(uniroot(equ,c(0,ti))$root)
}
t0Prediction2 <- function(X){
  X[1] -> ti; X[2] -> deltaV; X[3] -> aG; X[4] -> aM; X[5] -> bP; X[6] -> bM; X[7] -> theta; X[8] -> pdef
  equ = function(x) volume(x,aG,aM,bP,bM,theta,pdef) - (volume(ti,aG,aM,bP,bM,theta,pdef) - deltaV)
  return(uniroot(equ,c(0,ti))$root)
}
# aggregResults <- function(ti, om0, aG, aM, bP,bM,theta, Vext, )


#### Arguments ####
## tvec: maturity vector
## t1: first maturity
## om1: first (post-logging) proportion of timber species in the total volume
## omR: proportion of timber species in recruits' volume
## ag, am, bg, bm, th: volume model parameters
## int, slope: logging precision model parameters (intercept, slope)

omega_t = function(tvec, t1, om1, omR, ag, am, bg, bm, th, int, slope){
  
  t = seq.int(t1,max(tvec)); ord_t = floor(tvec - t1 + 1)
  
  V = ag/th*(1-(th*exp(-bg*t)-bg*exp(-th*t))/(th-bg))-am/th*(1-(th*exp(-bm*t)-bm*exp(-th*t))/(th-bm))
  
  dVG = ag*bg/(th-bg)*(exp(-bg*t)-exp(-th*t))+am*(1-(th*exp(-bm*t)-bm*exp(-th*t))/(th-bm))
  
  dVM = am*(1-exp(-bm*t))
  
  pR = 1/(1+exp(-(int+slope*log(V))))
  
  om = om1
  
  for (i in 1:(length(t)-1)) {
    om2 = min((om[i]*V[i] + dVG[i]*(pR[i]*omR + (1-pR[i])*om[i]) - dVM[i]*om[i])/V[i+1], 1)
    # if (V[i+1] == 0) om2 = omR
    om = c(om,om2)
  }
  
  return(om[ord_t])
}