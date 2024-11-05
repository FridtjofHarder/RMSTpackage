##########################################################################
###  Sample size calculation for an equivalence trial with binary endpoint
###     based on 
###     Makuch R, Simon R
###     Sample size requirements for evaluating a conservative therapy.
###     Cancer Treatment Reports [1978, 62(7):1037-1040]
###
###                                        DH 2015
#########################################################################

##  p0    true probability of the standard treatment 
##  p1    true probability of the reduced treatment treatment 
##  delta tolerance difference for (p0-p1)

MakuchSimon<-function(p0,p1,delta,alpha=.9,power=.8){
  Zalpha<-qnorm(alpha)
  Zbeta<-qnorm(power)
  
  n<-(p0*(1-p0)+p1*(1-p1))*((Zalpha+Zbeta)/(delta-(p0-p1)))^2
  print("For each treatment group")
  return(ceiling(n))
}

# MakuchSimon(.8,.8,.1,alpha=0.90) ## Gives same result as in Paper
# MakuchSimon(.8,.75,.1,alpha=0.90) ## Gives same result as in Paper



#---------------------------------------------------------------
MakuchSimon(.93,.93,0.055,alpha=.95,power=.80)*2/.9

