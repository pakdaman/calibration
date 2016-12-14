function [ out ] = predict( bbq, PTE, option )
% This function used for calibrating the probabilities
% Input: 
%       - bbq: the BBQ model 
%       - PTE : vector of Uncalibrated probabilities
%       - option: 0 use model selection, 1 use model Averaging   
% Output:
%       - out : vector of calibrated probabilities
    out = zeros(length(PTE),1);
    BBQModel = bbq.prunedModel; %it is also possible to use non-pruned version bbq.model
    if option == 1 % Use model Averaging
        parfor i=1:length(PTE)
            out(i,:) = getMA_local(BBQModel,PTE(i));
        end   
    elseif option == 0 % Use Model Selection
        out = getHistPr_local(BBQModel{1}.histModel, BBQModel{1}.cutPoints, PTE);
    end
end

function [res] = getMA_local(BBQModel, x)
    N = length(BBQModel);
    p = zeros(N,1);
    SV = BBQModel{1}.SV;  % It has already the relative likelihood
    for i=1:N
        p(i,:) = getHistPr_local(BBQModel{i}.histModel, BBQModel{i}.cutPoints, x);
    end
    
    res = SV'*p /sum(SV);    
end

function [res] = getHistPr_local(histModel, cutPoints, PTE)
    N = length(PTE);
    B = length(histModel);
    cutPoints = [0;cutPoints;1];
    res = zeros(N,1);
    
    for i=1:N
        x = PTE(i);
        minIdx = 1;
        maxIdx = B+1;
        
        while (maxIdx - minIdx)>1
            midIdx = floor((minIdx+maxIdx)/2);
            if x > cutPoints(midIdx)
                minIdx = midIdx;

            elseif x < cutPoints(midIdx)
                maxIdx = midIdx;
            else
                minIdx = midIdx;
                break;
            end

        end

        idx = minIdx;
        res(i,:) = histModel{idx}.P;    
        
        % For Handling Odd cases where there are more than one bin with exactly the
        % same min-max but different p (This is just happens for very rare cases in Naive Baye)

        cnt = 1;
        k = idx -1;
        while k>=1
            if (histModel{k}.min== histModel{idx}.min) && (histModel{k}.max==histModel{idx}.max)
                res(i,:) = res(i,:) + histModel{k}.P;
                k = k - 1;
                cnt = cnt +1;
            else
                break;
            end
        end
        
        k = idx +1;
        while k <= B
            if (histModel{k}.min== histModel{idx}.min) && (histModel{k}.max==histModel{idx}.max)
                res(i,:) = res(i,:) + histModel{k}.P;
                k = k + 1;
                cnt = cnt +1;
            else
                break;
            end
        end
        
        res(i,:) = res(i,:)/cnt;
        
    end
end