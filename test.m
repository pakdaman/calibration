clear; clc;
load data;

%% Apply linear the Logistic regression model to make predictions
disp('Using Linear SVM ')
svmStruct = svmtrain(XTR,YTR,'ShowPlot',true);
title('Linear SVM')
PTR = Mysvmclassify(svmStruct,XTR);
PTR = exp(PTR)./(1+exp(PTR));% Convert output of SVM to uncalibrated probs 
PTE = Mysvmclassify(svmStruct,XTE);
PTE = exp(PTE)./(1+exp(PTE));% Convert output of SVM to uncalibrated probs 

disp('Performance of (Linear) SVM probabilities : ')
M = getMeasures(PTE,YTE)
% Build BBQ Model
options.N0 = 2;
BBQ = build(PTR, YTR, options);
PTE_bbq = predict(BBQ, PTE, 1);
disp('Performance of Calibrated Probabilities using BBQ : ')
M_bbq = getMeasures(PTE_bbq,YTE)


%% Apply Quadratic the Logistic regression model to make predictions
disp('Using Quadratic SVM ')
figure;
svmStruct = svmtrain(XTR,YTR,'Kernel_Function','polynomial','Polyorder',2, 'ShowPlot',true);
title('Quadratic SVM')

PTR = Mysvmclassify(svmStruct,XTR);
PTR = exp(PTR)./(1+exp(PTR));% Convert output of SVM to uncalibrated probs 
PTE = Mysvmclassify(svmStruct,XTE);
PTE = exp(PTE)./(1+exp(PTE));% Convert output of SVM to uncalibrated probs 

disp('Performance of (Quadratic) SVM probabilities : ')
M = getMeasures(PTE,YTE)
% Build BBQ Model
options.N0 = 2;
BBQ = build(PTR, YTR, options);

% Calibrated the prediction using BBQ
PTE_bbq = predict(BBQ, PTE, 1);
disp('Performance of Calibrated Probabilities using BBQ : ')
M_bbq = getMeasures(PTE_bbq,YTE)

disp('End!')




