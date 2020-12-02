######################################################################
#  Population growth of structured population: Matrix Model 
#
#  Dec 1, 2020
#
######################################################################
#setwd("/YC_THU_MAC/Teaching_talk/Fall_2020/Biological_modeling/Lecture/Matrix_model/Lab")

## Create a Leslie matrix (based up Table 3.2, P62)

a=rep(0,16)

A=matrix(a, ncol=4)

A[1,]=c(1.6,1.5,0.25,0)

A[2,1]=0.8
A[3,2]=0.5
A[4,3]=0.25

## Checking A

## The initial age structure 

n0=c(50,50,50,50)
#n0=c(200,0,0,0)

## Time t1
n1=A%*%n0

lamda1=sum(n0)/sum(n1)

## Time t2
n2=A%*%n1
lamda2=sum(n2)/sum(n1)

### Population projections

nperm=10

pop.matrix=matrix(ncol=length(n0), nrow=nperm)

lamda=vector()

diff1=vector()
diff2=vector()
diff3=vector()
diff4=vector()

############ Starting point ####

nt=n0

for (i in 1:nperm)
{
  temp=nt
  nt=A%*%temp
  
  lamda[i]=sum(nt)/sum(temp)
  pop.matrix[i,]=nt
  
  diff1[i]=(nt[1]/sum(nt))-(temp[1]/sum(temp))
  diff2[i]=(nt[2]/sum(nt))-(temp[2]/sum(temp))
  diff3[i]=(nt[3]/sum(nt))-(temp[3]/sum(temp))
  diff4[i]=(nt[4]/sum(nt))-(temp[4]/sum(temp))
}

## Checking pop.matrix
pop.matrix

## Checking lamda 

lamda

## Calculate the eigenvalue

eigen_A=eigen(A)
eigen_A$values[1]

## Calculate the intrinsic rate of increase
r = log(eigen_A$values[1])

pdf(file="size_50.pdf")
{
  par(mfrow=c(1,1), pin=c(4,4))
  pop.lim=log10(max(pop.matrix))+0.3
  plot(seq(1, nperm), log10(pop.matrix[,1]), type="l", ylim=c(0,pop.lim), 
       xlab="Time", ylab="Log(N)", cex=3)
  lines(seq(1, nperm), log10(pop.matrix[,2]), col="blue", lty=2, cex=3)
  lines(seq(1, nperm), log10(pop.matrix[,3]), col="brown", lty=3, cex=3)
  lines(seq(1, nperm), log10(pop.matrix[,4]), col="darkgreen", lty=4, cex=3)
}
dev.off()


pdf(file="diff_50.pdf")
{
par(mfrow=c(2,2))
plot(seq(1, nperm), diff1, type="b", col="blue", ylim=c(-0.1,0.1), xlab="Time", ylab="Prop. Difference")
plot(seq(1, nperm), diff2, type="b", col="blue", ylim=c(-0.1,0.1), xlab="Time", ylab="Prop. Difference")
plot(seq(1, nperm), diff3, type="b", col="blue", ylim=c(-0.1,0.1), xlab="Time", ylab="Prop. Difference")
plot(seq(1, nperm), diff4, type="b", col="blue", ylim=c(-0.1,0.1), xlab="Time", ylab="Prop. Difference")
}
dev.off()

## Please repeate the process for n0=c(200,0,0,0) and compare the result

