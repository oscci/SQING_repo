#SQING_Repo script

#Power analysis for SQING
# by DVM Bishop 29th April 2017

require(pwr) # need to have package 'pwr' installed
# now specify Ns for 2 (recessive/dominant) or 3 (additive) genotype groups
# If only 2 groups, N3 set to zero
N1=100
N2=90
N3=0

#Specify value of power for which you want to compute effect size (r), 
#and value of r for which you want to compute power,  given this N

# For SQING we computed effect size (r) detectable with sample size, p = .05, with 80% power
myprepower<-.8
# As well as the power to detect an effect size (r) of .2 given the sample size
myprer<-.2

result1= pwr.r.test(n =N1+N2+N3 , sig.level =.05,power=myprepower,alternative="two.sided")
result2= pwr.r.test(n =N1+N2+N3 , sig.level =.05,r=myprer,alternative="two.sided")

paste("r detectable with this sample size, p = .05, and power of ",myprepower, " = ",result1$r)
paste("Power to detect r = ",myprer," at p = .05 with this sample size = ",result2$power)


