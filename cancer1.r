library(copent)
library(survival)
library(randomForestSRC)
library(fastcox)

# cancer data
data(cancer, package = "survival")
cancer = na.omit(cancer)
cancer$status = cancer$status - 1

# copula entropy
ce2 = 0
for(i in 4:10){
  ce2[i] = copent(cancer[,c(i,2)])
}

# survreg
surv2 = survreg(Surv(time,status)~., data = cancer)
pred2 = predict(surv2, cancer)
mae2a = mean(abs(pred2-cancer$time))
cidx2a = get.cindex(cancer$time,cancer$status,pred2)

# ce
survce2 = survreg(Surv(time,status)~sex+ph.ecog+ph.karno+pat.karno, data = cancer)
predce2 = predict(survce2, cancer)
mae2b = mean(abs(predce2-cancer$time))
cidx2b = get.cindex(cancer$time,cancer$status,predce2)

# random survival forest
rsf2 = rfsrc(Surv(time,status)~., data = cancer)
predrsf2 = predict(rsf2)
mae2c = mean(abs(predrsf2$predicted-cancer$time))
cidx2c = get.cindex(cancer$time,cancer$status,predrsf2$predicted)

# lasso-cox
ck2 = cocktail(as.matrix(cancer[,4:10]),cancer$time,cancer$status,nlambda = 30,alpha = 1)
predck2 = predict(ck2,as.matrix(cancer[,4:10]))
mae2d = mean(abs(predck2[,5]-cancer$time))
cidx2d = get.cindex(cancer$time,cancer$status,predck2[,5])

# plot variable importance
x11(width = 12, height = 12)
par(mfrow = c(2,2))
coef2 = surv2$coefficients[3:9]
barplot(coef2, ylab = "Coefficient value", xlab = "", main = "survreg")
bce2 = barplot(ce2[4:10], ylab = "Copula Entropy", ylim = c(min(ce2[4:10])-0.15,max(ce2[4:10])), main = "CE")
text(cex = 1, x = bce2, y = min(ce2[4:10])-0.1, labels = names(cancer)[4:10])
vs2 = var.select(rsf2,method = "md")
vs2a = vs2$varselect$depth[c(3,8,2,6,1,4,5)]
bce2 = barplot(vs2a, ylab = "Variable Importance", ylim = c(-0.3,max(vs2a)), main = "RSF")
text(cex = 1, x = bce2, y = -0.15, labels = names(cancer)[4:10])
ck2beta5 = ck2$beta[,5]
barplot(ck2beta5, ylab = "Coefficient value", main = "Lasso-Cox")

# plot prediction performance
x11(width = 10, height = 5)
par(mfrow = c(1,2))
mae2 = c(mae2a,mae2b,mae2c,mae2d)
names(mae2) = c("survreg","CE","RSF","Lasso-Cox")
bmae2 = barplot(mae2, ylab = "MAE", ylim = c(0, max(mae2)+20))
text(cex = 1, x = bmae2, y = mae2 + 12, labels = round(mae2,3), col = "red")
cidx2 = c(cidx2a,cidx2b,cidx2c,cidx2d)
names(cidx2) = c("survreg","CE","RSF","Lasso-Cox")
bcidx2 = barplot(cidx2, ylab = "C-Index", ylim = c(0, max(cidx2)+0.1))
text(cex = 1, x = bcidx2, y = cidx2+0.03, labels = round(cidx2,3), col = "red")
