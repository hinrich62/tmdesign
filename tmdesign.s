#R-code used to estimate the optimal allocation of tagged smolts between transported
#and in-river groups.
#inputs
#N -- the total number of tagged smolts
#SART -- SAR of transported juveniles
#SARM -- SAR of in-river juvenile migrants
tmdesign<-function(N=1000,SART=.01 ,SARM=.02,alpha=0.05){
if(N<=0){
warning("The number of juveniles released must exceed zero")
return(NULL)}
if((SART>1)|(SARM<0)){
warning("The SARs are survival rates and so must be less than one")
return(NULL)}
if((SARM>1)|(SARM<0)){
warning("The SARs are survival rates and so must be less than one")
return(NULL)}
s1<-sqrt((1-SART)/SART)
s2<-sqrt((1-SARM)/SARM)
ft<-s1/(s1+s2)
NT<-ft*N
NM<-N-NT
delta<-log(SART/SARM)
q<-qnorm(1-alpha/2)
tmvar<-(1-SART)/(SART*NT)+(1-SARM)/(SARM*NM)
se<-sqrt(tmvar)
power<-(1-pnorm(q*se,mean=delta,sd=se))+pnorm(-q*se,mean=delta,sd=se)
return(list(N=N,SART=SART, SARM=SARM ,alpha=alpha, NT=NT, NM=NM, delta=delta,se=se,cv=se/delta,power=power))
}
#outputs
#NT -- the optimal number of tagged transported juveniles
#NM – the optimal number of tagged in-river juvenile migrants
#delta is the true log(SART/SARM)
#se is the standard error of the estimate
#cv is the CV of the estimate
#power is the probability of rejecting the null hypothesis of delta=0.