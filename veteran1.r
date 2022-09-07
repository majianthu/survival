library(copent)
library(survival)
library(randomForestSRC)
library(fastcox)

data(cancer, package = "survival")
veteran$celltype = as.numeric(veteran$celltype)

ce3 = 0
for(i in c(1:2,5:8)){
  ce3[i] = copent(veteran[,c(i,3)])
}

# survreg
surv3 = survreg(Surv(time,status)~., data = veteran)
pred3 = predict(surv3, veteran)
mae3a = mean(abs(veteran$time-pred3))
cidx3a = get.cindex(veteran$time,veteran$status,pred3)

# ce
survce3 = survreg(Surv(time,status)~trt+celltype+karno+prior, data = veteran)
predce3 = predict(survce3, veteran)
mae3b = mean(abs(veteran$time-predce3))
cidx3b = get.cindex(veteran$time,veteran$status,predce3)

# random survival forest
rsf3 = rfsrc(Surv(time,status)~.,data = veteran)
predrsf3 = predict(rsf3)
mae3c = mean(abs(veteran$celltype-predrsf3$predicted))
cidx3c = get.cindex(veteran$time,veteran$status,predrsf3$predicted)

# regularized-cox
ck3 = cocktail(as.matrix(veteran[,c(1:2,5:8)]), veteran$time, veteran$status, nlambda = 20, alpha = 1)
predck3 = predict(ck3,as.matrix(veteran[,c(1:2,5:8)]))
mae3d = mean(abs(veteran$time-predck3[,4]))
cidx3d = get.cindex(veteran$time,veteran$status,predck3[,4])

# plot variable importance
x11(width = 12, height = 12)
par(mfrow = c(2,2))
coef3 = surv3$coefficients[2:7]
barplot(coef3, ylab = "Coefficient value", main = "survreg")
bce3 = barplot(ce3[c(1:2,5:8)], ylab = "Copula Entropy", ylim = c(min(ce3[c(1:2,5:8)])-0.2,max(ce3[c(1:2,5:8)])), main = "CE")
text(cex = 1, x = bce3, y = min(ce3[c(1:2,5:8)])-0.1, labels = names(veteran)[c(1:2,5:8)])
vs3 = var.select(rsf3,method = "md")
vs3a = vs3$varselect$depth[c(5,2,1,4,3,6)]
bce3 = barplot(vs3a, ylab = "Variable Importance", ylim = c(-0.3,max(vs3a)), main = "RSF")
text(cex = 1, x = bce3, y = -0.2, labels = names(veteran)[c(1:2,5:8)])
ck3beta5 = ck3$beta[,5]
barplot(ck3beta5, ylab = "Coefficient value", main = "Lasso-Cox")

# plot performane measures
x11(width = 10, height = 5)
par(mfrow = c(1,2))
mae3 = c(mae3a,mae3b,mae3c,mae3d)
names(mae3) = c("survreg","CE","RSF","Lasso-Cox")
bmae3 = barplot(mae3, ylab = "MAE", ylim = c(0, max(mae3)+10))
text(cex = 1, x = bmae3, y = mae3 + 5, labels = round(mae3,3), col = "red")
cidx3 = c(cidx3a,cidx3b,cidx3c,cidx3d)
names(cidx3) = c("survreg","CE","RSF","Lasso-Cox")
bcidx3 = barplot(cidx3, ylab = "C-Index", ylim = c(0, max(cidx3)+0.1))
text(cex = 1, x = bcidx3, y = cidx3 + 0.03, labels = round(cidx3,4), col = "red")
