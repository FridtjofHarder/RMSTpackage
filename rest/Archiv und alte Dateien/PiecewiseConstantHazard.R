

# Specification of piecwise constant hazards
CumPeriods <- c(3, 6, 12, 24, 30, 36)       # after this number of months
HazRates <- c(.01,.02, .05, .08,.04, .02)  

hazard<-function(t){ 
  if(t<=0 | t>max(CumPeriods)) res=0 else {
    res=rev(HazRates[(t < CumPeriods)])
    res=res[length(res)]}
  return(res)
  }
hazard<-Vectorize(hazard)
hazard(c(-1,3,13,25,35,47))

plot(hazard,xlim=c(-6,42))

Hazard<-function(t){integrate(hazard,lower=0,upper=t)$value}
# Better analytically: piecewise linear function
Hazard(0)

Hazard<-Vectorize(Hazard)
plot(Hazard,xlim=c(0,42))

fail=.7
S<-function(t) { exp(-Hazard(t))}

plot(S,xlim=c(0,42),xaxp=c(0,42,14),ylim=c(0,1))
abline(v=CumPeriods,lty=2)
# Is pretty flexible!


# ToDos:  
#    (1) Given CumPeriods (t1,..,tk) and 
#    Survival values S(t1)...S(tk) monoton decreasing.
#    Finde reasonably smooth interpolation by piecewise constant hazards
#
#    (2) Given S(t) given as piecewise constant hazards
#    Simulate efficently from this distribution. 
#
#    (3) Implement the formulas for RMST and RSDST in Royston 2013 Appendix
#        for S(t) given as piecewise constant hazards
#

