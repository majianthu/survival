library(copent)
library(randomForestSRC)
library(survival)
library(survsim)
library(fastcox)

# simulation
c6a = c7a = c8a = c9a = c10a = 0
c6c = c7c = c8c = c9c = c10c = 0
n1 = 1000 # number of subjects
sim.data <- simple.surv.sim(n=n1, foltime=100, 
                            dist.ev=c('weibull'), anc.ev=2, beta0.ev=1, 
                            anc.cens=0.85, beta0.cens=5,
                            # z=list(c("unif", 0.8, 1.2)), 
                            beta=list(c(1.4),c(1.2),c(0),c(1.2),c(0.2)), 
                            x=list(c("normal",0.4,1.1),c("normal",1,1.1),c("normal",0.7,1.1),c("normal",0.2,1.3),c("normal",0.2,1.1)))
# hist of time-to-event
x11(width = 10, height = 5)
hist(sim.data$stop,100, xlab = "days", main = "")

# random survival forest
rsf1 = rfsrc(Surv(stop,status)~x+x.1+x.2+x.3+x.4, data = sim.data)
vs1 = var.select(rsf1,method = "md")
vs1$varselect
vsrsf = vs1$varselect$depth[c(1,3,5,2,4)]
names(vsrsf) = c("x1","x2","x3","x4","x5")

# regularized-cox
ck1 = cocktail(as.matrix(sim.data[,6:10]),sim.data$stop,sim.data$status, nlambda = 20, alpha = 1)
betack = ck1$beta[,4]
names(betack) = c("x1","x2","x3","x4","x5")

# copula entropy
c6a = c6a + copent(sim.data[,c(4,6)])
c7a = c7a + copent(sim.data[,c(4,7)])
c8a = c8a + copent(sim.data[,c(4,8)])
c9a = c9a + copent(sim.data[,c(4,9)])
c10a = c10a + copent(sim.data[,c(4,10)])
# with censor mark
sim.data[,2] = sim.data[,2] + 0.00000001 * runif(n1)
c6c = c6c + copent(sim.data[,c(2,4,6)])
c7c = c7c + copent(sim.data[,c(2,4,7)])
c8c = c8c + copent(sim.data[,c(2,4,8)])
c9c = c9c + copent(sim.data[,c(2,4,9)])
c10c = c10c + copent(sim.data[,c(2,4,10)])

c1 = c(c6a,c7a,c8a,c9a,c10a)
names(c1) = c("x1","x2","x3","x4","x5")
c2 = c(c6c,c7c,c8c,c9c,c10c)
names(c2) = c("x1","x2","x3","x4","x5")

x11()
par(mfrow= c(2,2))
barplot(c1, main = "CE1")
barplot(c2, main = "CE2")
barplot(vsrsf, main = "RSF")
barplot(betack, main = "Lasso-cox")
