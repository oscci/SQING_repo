# Simulation of genotype-phenotype relations by DVM Bishop
# 22nd April 2017
# Updated by Paul Thompson to add error bars, August 2017
#-------------------------------------------------------------------------
setwd("~/Dropbox/Projects2016/SQING/SQING R scripts")
library(MASS) #for mvrnorm function to make multivariate normal distributed vars
library(yarrr)
library(tidyr)

options(scipen=999) #disable scientific notation.
#-------------------------------------------------------------------------
# Specify parameters to create correlated variables
#-------------------------------------------------------------------------
nVar<-2 #number of simulated variables 
myM<-0 # Mean score for simulated variables
myVar<-1 #Variance for simulated variables
myN<-300 #set sample size per group (You can vary this to see the effect)
i <-0
nSims<-100 #arbitary N simulations
mycutoffs<-c(-2.32,-.43,0,.43,.67)#range of cutoffs to try
mycentiles<-c(99,66,50,33,25)
ncutoff<-length(mycutoffs)
summarytable<-data.frame(matrix(rep(NA,11*ncutoff*nSims),ncol=11)) #initialise table to hold results
colnames(summarytable)<-c('Nsub','truer','cutoff','p.r','p.chi','N_aa','N_Aa','N_AA',
                          'mean_aa','mean_Aa','mean_AA')
for (k in 1:nSims){
  for (myCorr in c(.2)){ #correlation between variables; can loop through various values
    
    myCov<-matrix(rep(myCorr,nVar*nVar),nrow=nVar) #rep(x,y) means generate y values of x
    diag(myCov)<-rep(myVar,nVar) 
    #----------------------------------------------------------------------------------------
    # Generate a sample from a multivariate normal distribution with the specified correlation matrix
    #----------------------------------------------------------------------------------------
    mydata<-data.frame(mvrnorm(n = myN*100, rep(myM,nVar), myCov)) #make big population to select from
    colnames(mydata)<-c('myz','pheno')
    MAF<-.5 # minor allele frequency
    # convert the random normal deviates to genotype values of 0, 1 or 2, depending on MAF
    Naa<-MAF*MAF
    NAA<-(1-MAF)^2
    NaA<-1-Naa-NAA
    mydata$mygeno<-1
    critA<--qnorm(NAA)
    mydata[mydata$myz>critA,3]<-2
    crita<-qnorm(Naa)
    mydata[mydata$myz<crita,3]<-0
    
    #test significant of correlation of genotype/phenotype for full population
    rtest<-cor.test(mydata$pheno,mydata$mygeno)
    myr<-rtest$estimate
    
    regtest<-summary(lm(pheno~mygeno,data=mydata))
    myp<-regtest$coefficients[8] #pvalue for linear regression of pheno on geno
    myevalue<-c(as.integer(myN/4),as.integer(myN/2),as.integer(myN/4)) #expected distribution of genotypes
    myovalue<-table(mydata$mygeno[1:myN]) #observed distribution of genotypes
    mychi<-chisq.test(rbind(myevalue,myovalue))
    
    for (j in mycutoffs){ #range of cutoffs to try
      
      i<-i+1 #increase index for summary table row
      summarytable[i,1]<-myN
      summarytable[i,2]<-myr
      summarytable[i,3]<-j #cutoff on this run
      selsample <- mydata[mydata$pheno>j,]
      selsample<-selsample[1:myN,]
      rtest<-cor.test(selsample$pheno,selsample$mygeno)
      myp<-rtest$p.value #alternative way of getting p for regression
      if (myp< .001){myp <- .001} #set floor for low p-values to help plotting
      Nselected<-nrow(selsample)
      myevalue[2]<-as.integer(Nselected*NaA)
      myevalue[1]<-as.integer(Nselected*Naa)
      myevalue[3]<-as.integer(Nselected*NAA)
      myovalue<-table(selsample$mygeno)
      mychi<-chisq.test(rbind(myevalue,myovalue))
      myp2<-mychi$p.value
      if (myp2< .001){myp2<- .001} #set floor for low p-values to help plotting
      summarytable[i,4]<-myp
      summarytable[i,5]<-myp2
      summarytable[i,6:8]<-myovalue
      mymeans<-aggregate(selsample$pheno,by=list(selsample$mygeno),FUN=mean)
      summarytable[i,9:11]<-mymeans[,2]
    }
  }
}
#----------------------------------------------------------------------------------------
summarytable$log.pr<-(-log10(summarytable[,4]))
summarytable$log.pchi<-(-log10(summarytable[,5]))

myrange=1:length(mycutoffs)

#-------------------------------------------------------------
# aggregate data to get means by cutoff

aggdata<-aggregate(summarytable[,6:13], by=list(summarytable$cutoff),
                   FUN=mean, na.rm=TRUE)
colnames(aggdata)[1]<- 'cutoff'
aggdata$centiles<-mycentiles

#mysub<-paste('Population r =',myCorr,': Sample size =',myN)
# plot(aggdata$centiles,aggdata$log.pchi,ylim=c(0,2.5),xlim=rev(range(mycentiles)),main="",
#      xlab='% Ability-range included',ylab=expression(paste(-log[10], ' p-value')),xaxt='n',cex.axis=1.2,cex.lab=1.2)
# axis(side=1,at=mycentiles,cex.axis=1.2)
# lines(aggdata$centiles,aggdata$log.pchi)
# lines(aggdata$centiles,aggdata$log.pr,col='blue',type='o',pch=18)
# abline(h=-log10(.05),col='red',lty=2)
# text(53,.3,'(a) Chi-square test',cex=1.2) 
# text(60,1.8,'(b) Regression',cex=1.2) 
# text(90,1.15,'p = .05',cex=1.2)
# 
 colnames(aggdata)[5:7]<-c('aa','aA','AA')


#Illustrate how means and N change with cutoffs at 3 levels
# par(mfrow=c(1,3)) #set up plot frame with 3 plots in 1 row
# 
# barplot(as.matrix(aggdata[1,5:7]),main='No cutoff',ylim=c(-.2,1.25),ylab='zscore')
# text(.6,.04,paste('N = ',as.integer(aggdata[1,2])),cex=.8)
# text(1.8,.08,paste('N = ',as.integer(aggdata[1,3])),cex=.8)
# text(3.1,.35,paste('N = ',as.integer(aggdata[1,4])),cex=.8)
# 
# barplot(as.matrix(aggdata[3,5:7]),main='Cutoff = 0',ylim=c(-.2,1.25))
# text(.6,.8,paste('N = ',as.integer(aggdata[3,2])),cex=.8)
# text(1.8,.85,paste('N = ',as.integer(aggdata[3,3])),cex=.8)
# text(3.1,.95,paste('N = ',as.integer(aggdata[3,4])),cex=.8)
# 
# barplot(as.matrix(aggdata[4,5:7]),main='Cutoff = .5',ylim=c(-.2,1.25))
# text(.6,1.15,paste('N = ',as.integer(aggdata[4,2])),cex=.8)
# text(1.8,1.18,paste('N = ',as.integer(aggdata[4,3])),cex=.8)
# text(3.1,1.235,paste('N = ',as.integer(aggdata[4,4])),cex=.8)
# 


summarytable.cut1<-summarytable[summarytable$cutoff %in% c(-2.32),c(9:11)]

summarytable.cut2<-summarytable[summarytable$cutoff %in% c(0),c(9:11)]

summarytable.cut3<-summarytable[summarytable$cutoff %in% c(0.43),c(9:11)]

colnames(summarytable.cut1)[1:3]<-colnames(summarytable.cut2)[1:3]<-colnames(summarytable.cut3)[1:3]<-c('aa','aA','AA')

summarytable.cut2$cutoff<-rep("Top 50%",100)
summarytable.cut1$cutoff<-rep("No cutoff",100)
summarytable.cut3$cutoff<-rep("Top 33%",100)


plot.data<-rbind(summarytable.cut1,summarytable.cut2,summarytable.cut3)

long.plot.data<-gather(plot.data,genotype,z_score,aa:AA)

long.plot.data$cutoff<-factor(long.plot.data$cutoff,levels=c("No cutoff", "Top 50%", "Top 33%"))
long.plot.data$genotype<-factor(long.plot.data$genotype,levels=c("aa", "aA", "AA"))

test.pirate<-pirateplot(formula = z_score ~ genotype+cutoff,
           data = long.plot.data,
           xlab = "genotype",
           ylab = "z-score",
           main = "",
           theme = 2,  # Start with theme 2
           pal = "black",
           inf.f.o = 0, # Turn off inf fill
           inf.b.o = 0, # Turn off inf border
           point.o = 0,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           bar.b.o = 1,
           avg.line.o = 0, # Turn off average line
           point.col = "black",
           quant = c(.1, .9), # 10th and 90th quantiles
           quant.col = "black",plot=F) # Black points

tiff("Plot5.tiff", width = 6, height = 4, units = 'in', res = 300)
pirateplot(formula = z_score ~ genotype+cutoff,
           data = long.plot.data,
           xlab = "genotype",
           ylab = "z-score",
           main = "",
           theme = 2,  # Start with theme 2
           pal = "black",
           inf.f.o = 0, # Turn off inf fill
           inf.b.o = 0, # Turn off inf border
           point.o = 0,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           bar.b.o = 1,
           avg.line.o = 0, # Turn off average line
           point.col = "black",
           quant = c(.1, .9), # 10th and 90th quantiles
           quant.col = "black")


text(1,.1,'N = 75',cex=.6)
text(2,.3,'N = 150',cex=.6)
text(3,.5,'N = 75',cex=.6)

text(5,.9,'N = 60',cex=.6)
text(6,1,'N = 150',cex=.6)
text(7,1.1,'N = 90',cex=.6)

text(9,1.3,'N = 54',cex=.6)
text(10,1.3,'N = 149',cex=.6)
text(11,1.4,'N = 97',cex=.6)

dev.off()

#################################################
library(ggplot2)

long.agg<-gather(aggdata[,8:10],Test,logpval,log.pr:log.pchi)
long.agg$Test<-factor(long.agg$Test,levels=c("log.pr","log.pchi"))
levels(long.agg$Test)<-c("Regression","Chi-square test")

tiff("Plot6.tiff", width = 6, height = 4, units = 'in', res = 300)

ggplot(long.agg,aes(x=centiles,y=logpval,group=Test,color=Test,shape=Test))+geom_line(size=1)+geom_point(size=4)+ geom_hline(yintercept = -log10(.05),color="black",linetype="dashed",size=1) + annotate("text", label = "p = 0.05", x = 90, y = 1.2, size = 6)+xlab("% Ability-range included")+ylab(expression(paste(-log[10], ' p-value')))+theme_bw()+theme(legend.title=element_blank(),legend.position="top", legend.direction="horizontal",axis.text=element_text(size=14),axis.title=element_text(size=14),legend.text=element_text(size=14))+scale_x_continuous(trans = "reverse", breaks =c(25,33,50,66,99)) + scale_color_manual(values=c("darkslategray", "red"))+ scale_shape_manual(values = c(1, 19))
 dev.off()
par(mfrow = c(1,1))      
       
     xlab='% Ability-range included',ylab=expression(paste(-log[10], ' p-value')),xaxt='n',cex.axis=1.2,cex.lab=1.2)
axis(side=1,at=mycentiles,cex.axis=1.2)
lines(aggdata$centiles,aggdata$log.pchi)
lines(aggdata$centiles,aggdata$log.pr,col='blue',type='o',pch=18)
abline(h=-log10(.05),col='red',lty=2)
text(53,.3,'(a) Chi-square test',cex=1.2) 
text(60,1.8,'(b) Regression',cex=1.2) 
text(90,1.15,'p = .05',cex=1.2)