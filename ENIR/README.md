
Quick Start 
===========
This is a Matlab+R implementation of ENIR code. 
This is a temporary code and it will be replaced with a single R package enir soon.

How to use ENIR Matlab Code:

+ First Run test.m
% A sample code showing how to build a BBQ model and use it
% As a sanity check the output should match the output.txt

Two main Functions:
===========
+ build.m :
function [ enir ] = build( PTR, YTR, options )
% This function is used to build BBQ model
% Input:
%   PTR: Vector of predicted values
%   YTR: Vector of true labels {0,1}
%   Options: a Structure that includes optional choices
% Output:
%   enir: the ENIR moded

+ predict.m:
function [ out ] = predict( bbq, PTE, option )
% This function used for calibrating the probabilities
% Input: 
%       - bbq: the BBQ model 
%       - PTE : vector of Uncalibrated probabilities
%       - option: 0 use model selection, 1 use model Averaging   
% Output:
%       - out : vector of calibrated probabilities
