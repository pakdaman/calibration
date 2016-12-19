function [out,f] = Mysvmdecision(Xnew,svm_struct)
%SVMDECISION evaluates the SVM decision function

%   Copyright 2004-2006 The MathWorks, Inc.
%   $Revision: 1.1.12.4 $  $Date: 2006/06/16 20:07:18 $

sv = svm_struct.SupportVectors;
alphaHat = svm_struct.Alpha;
bias = svm_struct.Bias;
kfun = svm_struct.KernelFunction;
kfunargs = svm_struct.KernelFunctionArgs;

f = (feval(kfun,sv,Xnew,kfunargs{:})'*alphaHat(:)) + bias;
% out = sign(f);
% % points on the boundary are assigned to class 1
% out(out==0) = 1;

out = f;