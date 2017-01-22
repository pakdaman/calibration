# install.packages("e1071", dependencies = TRUE)
# install.packages("R.matlab", dependencies = TRUE)

library(R.matlab)
library(e1071)
library(ELiTE)

data = readMat("data.mat")
XTrain = data$XTR
XTest = data$XTE

ztr = data$YTR
zte = data$YTE

# Plot training data
plot(XTrain,col=ifelse(ztr>0,1,2))
legend("topleft",c('Positive','Negative'),col=seq(2),pch=1,text.col=seq(2))

svm.model = e1071::svm(XTrain, ztr, decision.values = TRUE, kernel = 'linear', cost=1)
summary(svm.model)
ytr = stats::predict(svm.model, XTrain)
ytr = exp(ytr)/(1+exp(ytr));
elite.model = elite.build(ytr, ztr, 'AICc')
yte = stats::predict(svm.model, XTest)
yte = exp(yte)/(1+exp(yte));
yte.cal =  elite.predict(elite.model, yte, 1);

print("Evaluation measures for Linear SVM:")
Mobj = elite.getMeasures(yte,zte)
print(Mobj)
print("Evaluation measures for Calibrated Linear SVM using ENIR:")
Mobj.cal = elite.getMeasures(yte.cal,zte)
print(Mobj.cal)



svm.model = e1071::svm(XTrain, ztr, decision.values = TRUE, kernel = 'polynomial', d = 2)
summary(svm.model)
ytr = stats::predict(svm.model, XTrain)
ytr = exp(ytr)/(1+exp(ytr));
elite.model = elite.build(ytr, ztr, 'AICc')
yte = stats::predict(svm.model, XTest)
yte = exp(yte)/(1+exp(yte));
yte.cal =  elite.predict(elite.model, yte, 1);

print("Evaluation measures for Quadratic SVM:")
Mobj = elite.getMeasures(yte,zte)
print(Mobj)
print("Evaluation measures for Calibrated Quadratic SVM using ENIR:")
Mobj.cal = elite.getMeasures(yte.cal,zte)
print(Mobj.cal)

