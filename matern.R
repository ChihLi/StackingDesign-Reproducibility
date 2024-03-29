### matern kernel function ###
matern.kernel <- function(r,k,derivative=0){
  # r is the distance
  # k is the smoothness parameter
  
  if(k==1/2){ #nu=0.5
    if(derivative == 0) out <- exp(-r)
    if(derivative == 1) out <- -exp(-r)
  }else if(k==3/2){ #nu=1.5
    if(derivative == 0) out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
    if(derivative == 1) out <- -3*r*exp(-sqrt(3)*r)
  }else if(k==5/2){ #nu=2.5
    if(derivative == 0) out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
    if(derivative == 1) out <- -(r*(5^(3/2)*r+5)*exp(-sqrt(5)*r))/3
  }else if(k==7/2){ #nu=3.5
    if(derivative == 0) out <- (1+r*sqrt(7)+2.8*r^2+7/15*sqrt(7)*r^3) * exp(-r*sqrt(7))
    if(derivative == 1) out <- -(r*(49*r^2+3*7^(3/2)*r+21)*exp(-sqrt(7)*r))/15
  }else if(k==9/2){ #nu=4.5
    if(derivative == 0) out <- (1+r*sqrt(9)+27*r^2/7+18/7*r^3+27/35*r^4) * exp(-r*3)
    if(derivative == 1) out <- -((81*r^4+162*r^3+135*r^2+45*r)*exp(-3*r))/35
  }
  
  return(out)
}