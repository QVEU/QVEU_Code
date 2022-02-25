#powerAnalysis
#Author: Patrick Dolan
#Simple worksheet for dealing with Powe

SD=1.8
a=0.1
b=0.05
power=1-b

effect_sizes = c(0.3, 0.5, 0.7)
ChangeValues = effect_sizes*SD
sample_sizes = seq(3,20)

Tt<-function(effect,sample) {
  t=sqrt(sample)*effect
  Pval=pt(t,df = 2*(sample)-2)
  return(Pval)
}

plot(type="l",sample_sizes,Tt(.7,sample_sizes), ylab="Power",xlab = "Sample Size",ylim=c(0.2,1),col="blue",lwd=2)
lines(sample_sizes,Tt(.5,sample_sizes),lwd=2)
lines(sample_sizes,Tt(.3,sample_sizes),col="darkred",lwd=2)
text(x = 15,y=0.75,"Power: 0.8",col="red")
abline(h = 0.8, col="red", lty=2)
legend(12,0.4,c("0.3 (0.54)","0.5 (0.9)","0.7 (1.26)"),lty = 1,col=c("darkred","black","blue"),cex = 0.6,title = "Eff. Size (Clin. Sc.)")

