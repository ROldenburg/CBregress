# CBregression R
# Reinhard Oldenburg - roldenburg@gmx.de 
# Linear Regression for the model y=a*x+b

library("nloptr")

# orthogonal regression - returns a vector of the slope and intercept (a,b)
orthoreg<-function(x,y){
  a=(var(y)-var(x)+sqrt(4*cov(x,y)^2+(var(y)-var(x))^2))/(2*cov(x,y))
  return(c(a,mean(y)-a*mean(x)) ) }

# case based regression
CBregress=function(x,y) {
  n=length(x); if(length(y)!=n) stop("vectors of unequal length")
  mx=mean(x); my=mean(y); x=x-mx; y=y-my;
  result=list()
  p1=1000;  p2=1000; p3=1000; p4=10; # penalty weights
  cxy=cov(x,y); vx=var(x); vy=var(y)
  # var v = c(L1,...,Ln,a)
  F=function(v){
    L=v[1:n]
    return(sum((x-L)^2) +  sum((y-L*v[n+1])^2) 
           +p1*sum((y-v[n+1]*L)*(x-L) )^2/max(0.00001,sqrt(sum((x-L)^2)*sum((y-L*v[n+1])^2)))
           +p2*sum((y-v[n+1]*L)*L)^2/max(0.00001,sqrt(sum((y-L*v[n+1])^2))*sqrt(sum(L^2)))+
             +p3*sum(L*(x-L))^2/max(0.00001, sqrt(sum((x-L)^2))*sqrt(sum(L^2))) 
           +p4*sum(L)^2 
    )    }
  H=function(v){return(-cxy*v[n+1])}
  v=c((x+y/cxy*vx)/2, cxy/vx)
  opts <- list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-9,maxeval=5000)
  res <- nloptr( x0=v, eval_f=F, eval_g_ineq=H, opts=opts) 
  vs=res$solution; 
  X0=vs[1:n]; result$X0=X0; result$a=a=vs[n+1]; result$b=my-a*mx
  result$v0=var(X0); result$vd=var(x-X0); result$ve=var(y-X0*a);
  return(result)
}


# Testcode : Performs a comparison between standard, orthogobal and case based regression
testRegs=function(m,n,A,sig1,sig2,generator=runif) {
  # m number of simulations
  # n sample size
  # A true regression coefficient
  # generator function like runify such that generator(n) gives random vector
  # calculates a vector of std. reg coefficient, orthogonal coefficient, case based coefficient
  # for each simulation, prints standard deviation of errors (i.e. estimates-true A), return variances 
  AW=c(0,0,0)
  for(k in 1:m) {
    x=generator(n); y=A*x+rnorm(n,0,sig1); x=x+rnorm(n,0,sig2)
    res=CBregress(x,y)
    aW=(c(cov(x,y)/var(x),orthoreg(x,y)[1],res$a)-A)^2
    AW=AW+aW
  }
  print(sqrt(AW/m))
  return(AW/m)
}

