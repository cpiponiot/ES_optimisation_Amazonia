f_cemi = function(t, vext=15, damage=15, acs=200, flwn = 0.87, p1=0.1,l1=0.1,l2=0.2){
  
  acs_left = damage + 0.12*acs
  decay_lwn = flwn*(t -  p1/l1*(1-exp(-l1*t)) - (1-p1)/l2*(1-exp(-l2*t)))
  decay_fwn = (1-flwn)*(t - (1-exp(-0.19*t))/0.19)
  
  acs_ext =  0.725*vext
  decay_sawmill = 0.67*t + 0.33*(t-(1-exp(-log(2)/30*t))/(log(2)/30))
  
  acs_left*(decay_lwn + decay_fwn) + acs_ext*decay_sawmill
}


f_recov = function(t, loss = 25, acs0 = 200, aSg = 120, aRr = 0.3, aRg = 1.5, 
                   bSg = 0.02, bSm = 0.007, bRr = 0.2, bRg = 0.05, bRm = 0.008, 
                   eta = 2, integral = FALSE){
  
  if (integral){
    cSg <- aSg*(t - (1-exp(-t*bSg))/bSg)
    cSm <- (aSg + acs0*(1-loss/100))*(t - (1-exp(-t*bSm))/bSm)
    cRr <- aRr*(t^2/2 + eta * (t-(1-exp(-t*bRr))/bRr)/bRr)
    cRg <- aRg*(t^2/2 - (t-(1-exp(-t*bRg))/bRg)/bRg)
    cRm <- (aRr+aRg)*(t^2/2 - (t-(1-exp(-t*bRm))/bRm)/bRm)
    
    return(cSg - cSm + cRr + cRg - cRm)
    
  } else {
    
    cSg <- aSg*(1-exp(-t*bSg))
    cSm <- (aSg + acs0*(1-loss/100))*(1-exp(-t*bSm))
    cRr <- aRr*(t + eta* (1-exp(-t*bRr))/bRr)
    cRg <- aRg*(t - (1-exp(-t*bRg))/bRg)
    cRm <- (aRr+aRg)*(t - (1-exp(-t*bRm))/bRm)
    
    return(cSg - cSm + cRr + cRg - cRm)
    
  }
  
}