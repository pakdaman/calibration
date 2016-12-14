Preface 
===========
This is a Matlab+R implementation of ENIR code. This is a research code and there is not guarantee on using the code.
This code will be replaced with a single R package "enir" soon.
To run the code you should first install the neariso R package referenced inside our icdm paper. 
For the sake of completeness, a version of neariso R package that is used in our experiments is also provided in this repository (neariso_1.0.tar.gz).
The code runs properly on R version 3.1.2 (2014-10-31) -- "Pumpkin Helmet"
 and Matlab version 8.4.0.150421 (R2014b)


Quick Start 
===========
How to use ENIR Matlab Code:
+ First Run test.m
% A sample code showing how to build a ENIR model and use it


Two main Functions:
===========
+ build.m :
function [ enir ] = build( PTR, YTR, options )
% This function is used to build ENIR model
% Input:
%   PTR: Vector of predicted values
%   YTR: Vector of true labels {0,1}
%   Options: a Structure that includes optional choices
% Output:
%   enir: the ENIR moded

+ predict.m:
function [ out ] = predict( enir, PTE, option )
% This function used for calibrating the probabilities
% Input: 
%       - enir: the ENIR model 
%       - PTE : vector of Uncalibrated probabilities
%       - option: 0 use model selection, 1 use model Averaging   
% Output:
%       - out : vector of calibrated probabilities

References 
===========

Mahdi Pakdaman Naeini, Gregory F. Cooper. "Binary Classifier Calibration using an Ensemble of Near Isotonic Regression Models”, 
IEEE Internation Conference on Data Mining (ICDM), Barcelona, Spain, December 2016.

