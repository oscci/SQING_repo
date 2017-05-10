#scratchpad

mytempx<-matrix(rep(NA,2*71),ncol=2)
mytempx[1:70,1]<-colnames(allstudy)# for checking col number of relevant columns
neworder<-c(3,4,54,13,5,6,7,8,9,29,30,10,33,34,39,11,25,27,50,12,16,44,22,15,45,38,41,28,31,35,37,17,18,48) #Authors omitted
mytempl<-matrix(rep(NA,2*34),ncol=2)

mytempl[1:34,1]<-c('Title', 'DOI', 'Journal','Conclusion','Sample', 'Discovery Sample N',
                   'Replication Sample N',
                   'Genes', 'Polymorphisms', 'N Polymorphisms','Correlated Polymorphisms',
                   'Phenotypes', 'N Behav. Phenotypes','Correlated Behav. Phenotypes',
                   'N Neuro Phenotypes','Analysis Method', 'Genetic models',
                   'N genetic models', 'Neuro methods','Selected Result',
                   'Selected Result Source','Ns for Selected Result',
                   'Comment on Selected Result','Selected Result Effect Size',
                   'Quasi Effect Size?','Imaging ROI approach',
                   'Mention of multiple comparisons','Correction for N genetic models',
                   'Correction for N polymorphisms','Correction for N behavioural phenotypes',
                   'Correction for N neuro phenotypes','Effect Size (r) Detectable with 80% Power',
                   'Power to Detect Effect Size (r) of .1','Author response to contact')
i=1
thisstudy=paste('Study',i)
#print(thisstudy)
colnames(mytempl)=c('Information',thisstudy)
mytempl[,2]<-t(allstudy[i,neworder])



temp=rep(0,30) #create vector to denote quasi effect sizes
myindex<-c(which(allstudy$Quasi_ES==0),which(allstudy$Quasi_ES=='n/a'))
mysubset<-allstudy[myindex,]
temp[myindex]<-1
plot(allstudy$SelRes_N,allstudy$SelRes_ES, main="Log sample size by Effect size (r)",log='x',
     xlab="Sample size ", ylab="Effect size (r)",
     pch=1+(15*temp))
