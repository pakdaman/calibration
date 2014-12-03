
Quick Start 
===========
How to use BBQ Matlab Code:

+ Run test.m
% A sample code showing how to build a BBQ model and use it
% As a sanity chech the output should match the output.txt

Two main Functions:
===========
+ build.m :
function [ bbq ] = build( PTR, YTR, options )
% This function is used to build BBQ model
% Input:
%   PTR: Vector of predicted values
%   YTR: Vector of true labels {0,1}
%   Options: a Structure that includes optional choices
% Output:
%   bbq: the BBQ moded

+ predict.m:
function [ out ] = predict( bbq, PTE, option )
% This function used for calibrating the probabilities
% Input: 
%       - bbq: the BBQ model 
%       - PTE : vector of Uncalibrated probabilities
%       - option: 0 use model selection, 1 use model Averaging   
% Output:
%       - out : vector of calibrated probabilities
