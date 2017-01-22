Preface 
===========
This is a research code and there is not guarantee on using the code.
This code includes an R package "ELiTE" that implements our ELiTE calibration method in our SDM 2016 paper.
To run the code you should first install the glmgen  package referenced inside our sdm paper. 
For the convenience, a version of glmgen R package that is used in our experiments is 
also provided in this repository (glmgen_elite.tar.gz). You should then install the ELiTE package which is also included in this
repository (ELiTE_1.0.tar.gz).

You can then run the sample.R file for an example code of using the enir package.
In order to run the sample.R, you should install the "R.matlab", and "e1071" in R using the following commands:

install.packages("e1071", dependencies = TRUE)
install.packages("R.matlab", dependencies = TRUE)

Two main methods + a utility function:
===========
+ elite.build(yTrain, zTrain);

--Input:
yTrain: vector of uncalibrated classification scores that are real numbers in the interval of [0,1]

zTrain: vector of corresponding true class. 1 indicates positive class and 0 indicates negative class.

option: 'BIC' (default) , 'AIC', or 'AICc' scoring functions (In the paper we have used AICc).

--Output:

A list of parameters that indicate the elite model



+ elite.predict(enirModel, yTest, option);

--Input: 

eliteModel: a ELiTE calibration model generated using the "elite.build()" method

yTest: vector of uncalibratd classification scores

option: set it to 1 (default) for running model averaging, and to 0 for running model selection

--output:

res: corresponding calibrated scores obtained using the elite calibration model


+ elite.getMeasures(y, z);

--Input:

y: vector of predictions (classification scores) which is in the interval [0, 1]

z: vector of true class of instances {0,1}

--Output: 

a list of evaluation measures including (RMSE, AUC, ACC, MCE, ECE)

References 
===========

Mahdi Pakdaman Naeini, Gregory F. Cooper. "Binary Classifier Calibration using an Ensemble of Linear Trend Estimation”, 
SIAM Data Mining (SDM), Miami, FL, 2016.

