function [res] = getMeasures(PTE, YTE )
%GETMEASURES Summary of this function goes here
%   Detailed explanation goes here
    if sum(isnan(PTE))>0
        disp('there are some nan value in predictions')
    end
    
    if sum(isnan(YTE))>0
        disp('there are some nan value in the labels')
    end
    
    idx = isnan(YTE)|isnan(PTE);
    YTE = YTE(~idx,:); PTE = PTE(~idx,:);

    res.RMSE = getRMSE(YTE,PTE);
    res.AUC = getAUC(YTE,PTE);
    res.ACC = 1-sum(YTE~=(PTE>=0.5))/length(YTE);    
    res.MCE = getMCE(YTE,PTE); % Computing the Max Calibratio Error among all bins 
    res.ECE = getECE(YTE,PTE); % Computing Average Calinration Error on different binns
end

