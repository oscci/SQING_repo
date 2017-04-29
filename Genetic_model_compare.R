#SQING_Repo script

# Genetic model correction simulation
# Testing how far recessive and dominant analyses are independent from additive
# And how inclusion of all 3 models affects false positive rate

# Simulation by DVM Bishop
# 5th April 2017



require(MASS) #for multivariate normal generation
require(data.table)
options(scipen=999) #turn off sci notation
options(digits=4)

#--------------------------------------------------------
# create dummy dataset
#--------------------------------------------------------
ndata=20000 #N simulated rows of data (subjects)
mydata <- data.frame(matrix(0, nrow = ndata, ncol = 6))
colnames(mydata)=c('allele1','allele2','Nmaj','pheno','recmodel','dommodel')
nrun=1000
myp=data.frame(matrix(0,nrow=nrun,ncol=4)) #cols for p values for each model, and last one for constellation
for (thisrun in 1:nrun){
  for (i in 1:2){
    mydata[,i]=runif(ndata) #random uniform number 0-1; used to select probabilistically
  }
  mydata[,4]=rnorm(ndata) #random normal deviate: used to simulate phenotype
  
  MAF=.5 #Minor allele frequency
  
  # Now convert probabilities in cols 1-2 to 1 and 0
  realeffect=0.0 #Set to zero if no true effect of gene on phenotype
  mydata$allele1[mydata$allele1< MAF] <- 0
  mydata$allele1[mydata$allele1>=MAF] <- 1
  mydata$allele2[mydata$allele2< MAF] <- 0
  mydata$allele2[mydata$allele2>=MAF] <- 1
  mydata$Nmaj=mydata$allele1+mydata$allele2
  
  # Now allocate columns 6 and 7 to 1 and 0 to capture contrast for specific model
  mydata$recmodel=mydata$Nmaj
  mydata$dommodel=mydata$Nmaj-1
  mydata$recmodel[mydata$recmodel==2]=1
  mydata$dommodel[mydata$dommodel<0]=0
  mydata$pheno=mydata$pheno+realeffect*mydata$Nmaj
  
  criticalp<-.05
  
  # Fit additive model, where phenotype increases linearly with N A alleles
  addfit <-lm(pheno ~ Nmaj, data=mydata)
  x=summary(addfit)
  myp[thisrun,1]=x$coefficients[2,4] #extract relevant p-values
  if (x$coefficients[2,4]<criticalp) {myp[thisrun,4]=myp[thisrun,4]+100}
  
  # Fit recessive model, where phenotype differs for 2 vs 0-1 alleles
  recfit <-lm(pheno ~ recmodel, data=mydata)
  x=summary(recfit)
  myp[thisrun,2]=x$coefficients[2,4]
  if (x$coefficients[2,4]<criticalp) {myp[thisrun,4]=myp[thisrun,4]+10}
  
  # Fit dominant model, where phenotype differs for 1-2 vs 0 alleles
  domfit <-lm(pheno ~ dommodel, data=mydata)
  x=summary(domfit)
  myp[thisrun,3]=x$coefficients[2,4]
  if (x$coefficients[2,4]<criticalp) {myp[thisrun,4]=myp[thisrun,4]+1}
  
}
colnames(myp)<-c('add_p','rec_p','dom_p','sigp')
#final column gives pattern of significant findings

allz=NA
allp=which(myp[,4]>0) 

overallp=length(allp)/nrun #how many runs give a significant value

my1=length(which(myp[,4]==1)) #only dominant model significant
my10=length(which(myp[,4]==10)) #only recessive model significant
my11=length(which(myp[,4]==11)) #only dom and recessive model significant
my100=length(which(myp[,4]==100)) #only additive model significant
my101=length(which(myp[,4]==101)) #additive and dominant significant
my110=length(which(myp[,4]==110)) #additive and recessive significant
my111=length(which(myp[,4]==111))  #all models significant
paste('True effect of gene:',realeffect)
paste('Dominant only significant',my1/nrun)
paste('Recessive only significant',my10/nrun)
paste('Additive only significant',my100/nrun)
paste('Recessive and dominant significant', my11/nrun)
paste('Additive and dominant significant',my101/nrun )
paste('Additive and recessive significant',my110/nrun)
paste('All three models significant',my111/nrun)
paste('At least one model significant',overallp)

# This simulation shows that to achieve overall 5% false positive rate when
# all 3 models are included, criticalp needs to be set to .025, rather than .05

# Plot a pie chart
par(mar=c(.2,.2,.2,.2))
x <- c(my111,my1,my10,my100,my101,my110)
labels <- c('All models','Dominant only','Recessive only','Additive only','Additive+Dominant','Additive\n+Recessive')

# Give the chart file a name (for plot in an external file)
png(file = "modelpie.jpg")

# Plot the chart.
pie(x,labels,cex=1.5)

# Save the file.
dev.off() #only needed if plot sent to external file