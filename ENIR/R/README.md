Preface 
===========
This is a research code and there is not guarantee on using the code.
This code includes an R package "enir" that implements our ENIR calibration method in our ICDM 2016 paper.
To run the code you should first install the neariso R package referenced inside our icdm paper. 
For the sake of completeness, a version of neariso R package that is used in our experiments is 
also provided in this repository (neariso_1.0.tar.gz). You should then install the enir package which is also included in this
repository (enir_1.0.tar.gz).

You can then run the sample.R file for an example code of using the enir package.
In order to run the sample.R, you should install the "R.matlab", and "e1071" in R using the following commands:

install.packages("e1071", dependencies = TRUE)
install.packages("R.matlab", dependencies = TRUE)

Two main methods + a utility function:
===========
+ enir.build(yTrain, zTrain);

--Input:
yTrain: vector of uncalibrated classification scores that are real numbers in the interval of [0,1]

zTrain: vector of corresponding true class. 1 indicates positive class and 0 indicates negative class.

option: 'BIC' (default) , 'AIC', or 'AICc' scoring functions.

--Output:

A list of parameters that indicate the enir model



+ enir.predict(enirModel, yTest);

--Input: 

enirModel: a ENIR calibration model generated using the "enir::build()" method

yTest: vector of uncalibratd classification scores


--output:

res: corresponding calibrated scores obtained using the ENIR calibration model


+ enir.getMeasures(y, z);

--Input:

y: vector of predictions (classification scores) which is in the interval [0, 1]

z: vector of true class of instances {0,1}

--Output: 

a list of evaluation measures including (RMSE, AUC, ACC, MCE, ECE)

References 
===========

Mahdi Pakdaman Naeini, Gregory F. Cooper. "Binary Classifier Calibration using an Ensemble of Near Isotonic Regression Models”, 
IEEE Internation Conference on Data Mining (ICDM), Barcelona, Spain, December 2016.

