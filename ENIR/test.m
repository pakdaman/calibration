clear; clc;
load data;

%% Apply linear SVM regression model to make predictions
disp('Using Linear SVM ')
svmStruct = svmtrain(XTR,YTR,'ShowPlot',true);
title('Linear SVM')
PTR = Mysvmclassify(svmStruct,XTR);
PTR = exp(PTR)./(1+exp(PTR));% Convert output of SVM to uncalibrated probs 
PTE = Mysvmclassify(svmStruct,XTE);
PTE = exp(PTE)./(1+exp(PTE));% Convert output of SVM to uncalibrated probs 

disp('Performance of (Linear) SVM probabilities : ')
M = getMeasures(PTE,YTE)
% Build ENIR Model
ENIR = build(PTR, YTR);
PTE_enir = predict(ENIR, PTE, 1);
disp('Performance of Calibrated Probabilities using BBQ : ')
M_enir = getMeasures(PTE_enir,YTE)


%% Apply Quadratic SVM regression model to make predictions
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
ENIR = build(PTR, YTR);

% Calibrated the prediction using BBQ
PTE_enir = predict(ENIR, PTE, 1);
disp('Performance of Calibrated Probabilities using BBQ : ')
M_enir = getMeasures(PTE_enir,YTE)

disp('End!')




