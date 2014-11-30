function [ res ] = getRMSE( Y, P )
%RMSE Summary of this function goes here
%   Detailed explanation goes here
    res = (((Y-P)'*(Y-P))/length(Y))^0.5;

end

