# Simulation of genotype-phenotype relations by DVM Bishop
# 22nd April 2017
#-------------------------------------------------------------------------
setwd("~/Dropbox/Projects2016/SQING/SQING R scripts")
library(MASS) #for mvrnorm function to make multivariate normal distributed vars

options(scipen=999) #disable scientific notation.
#-------------------------------------------------------------------------
# Specify parameters to create correlated variables
#-------------------------------------------------------------------------
nVar<-2 #number of simulated variables 
myM<-0 # Mean score for simulated variables
myVar<-1 #Variance for simulated variables
myN<-300 #set sample size per group (You can vary this to see the effect)
i <-0
nSims<-100 #arbitary N simulations (set to 100 for testing, 10000 for final run)
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
summarytable$log.pr<-log10(summarytable[,4])
 summarytable$log.pchi<-log10(summarytable[,5])

 myrange=1:length(mycutoffs)

#-------------------------------------------------------------
 # aggregate data to get means by cutoff 
 aggdata<-aggregate(summarytable[,6:13], by=list(summarytable$cutoff),
                    FUN=mean, na.rm=TRUE)
 colnames(aggdata)[1]<- 'cutoff'
 aggdata$centiles<-mycentiles

 #plot p-values for chi square/regression
 par(mar=c(1,1,0,0)) #set margins
 #set up file to write figure to (won't see it on screen)
 myfolder<-"~/Dropbox/Projects2016/SQING/draft SQING article/Figures/"
 myfilename<-'Figure6.png'
 png(paste(myfolder,myfilename,sep=''),
     width = 800, height = 620, units = "px", pointsize = 28)
 
#myhead<-'Log p-values for \n regression of phenotype on genotype'
 #don't include title for figure for article
#mysub<-paste('Population r =',myCorr,': Sample size =',myN)
 #this is subheading that would appear below x-axis
plot(aggdata$centiles,aggdata$log.pchi,ylim=c(-2.5,0),xlim=rev(range(mycentiles)),
     xlab='% Ability-range included',ylab='log10 p-value',xaxt='n')
axis(side=1,at=mycentiles)
lines(aggdata$centiles,aggdata$log.pchi)
lines(aggdata$centiles,aggdata$log.pr,col='blue',type='o',pch=16)
abline(h=log10(.05),col='red',lty=2)
text(53,-.3,'(a) Chi-square test') 
text(60,-2,'(b) Regression') 
text(90,-1.25,'p = .05')
colnames(aggdata)[5:7]<-c('aa','aA','AA')
dev.off() #turn off external plotting

#-----------------------------------------------------------------------
#Illustrate how means and N change with cutoffs at 3 levels
#-----------------------------------------------------------------------
myfilename<-'Figure5.png'
png(paste(myfolder,myfilename,sep=''),
    width = 960, height = 420, units = "px", pointsize = 28)

par(mfrow=c(1,3)) #set up plot frame with 3 plots in 1 row

barplot(as.matrix(aggdata[1,5:7]),main='No cutoff',ylim=c(-.2,1.35),
        ylab='z-score',col='yellow')
text(.65,.09,as.integer(aggdata[1,2]),cex=.85)
text(1.8,(aggdata[1,6]+.11),as.integer(aggdata[1,3]),cex=.85)
text(3.1,(aggdata[1,7]+.11),as.integer(aggdata[1,4]),cex=.85)
#NB text positions for N labels done by trial and error
# for +ve values, add around .11 to value

barplot(as.matrix(aggdata[3,5:7]),main='Top 50%',ylim=c(-.2,1.35),col='orange')
text(.65,(aggdata[3,5]+.11),as.integer(aggdata[3,2]),cex=.85)
text(1.8,(aggdata[3,6]+.11),as.integer(aggdata[3,3]),cex=.85)
text(3.1,(aggdata[3,7]+.11),as.integer(aggdata[3,4]),cex=.85)

barplot(as.matrix(aggdata[4,5:7]),main='Top 33%',ylim=c(-.2,1.35),col='red')
text(.65,(aggdata[4,5]+.11),as.integer(aggdata[4,2]),cex=.85)
text(1.8,(aggdata[4,6]+.11),as.integer(aggdata[4,3]),cex=.85)
text(3.1,(aggdata[4,7]+.11),as.integer(aggdata[4,4]),cex=.85)
dev.off() #turn off external plotting
