% Change for to parfor
function [ pout ] = predict( enir, pin, option )

% This function used for calibrating the probabilitiesusing ENIR model
% Input: 
%       - ENIR: the ENIR model learnt by lnEFBB
%       - pin : vector of Uncalibrated probabilities
%       - option: 0 use model selection, 1 use model Averaging   
% Output:
%       - out : vector of calibrated probabilities
%     pout = zeros(length(pin),1);
%     enirModel = enir.prunedModel;
    enirModel = enir.model;
    pout = zeros(length(pin),1);
    if option == 1 % Use model Averaging
        parfor i=1:length(pin)
            pout(i,:) = getMA_local(enirModel,pin(i));
        end   
    elseif option == 0 % Use Model Selection
        nirModel.y = enirModel.y;
        nirModel.uy = enirModel.uy;
        nirModel.z = enirModel.z;
        nirModel.df = enirModel.df(1);
        nirModel.beta = enirModel.beta(:,1);
        nirModel.lambda = enirModel.lambda(1);
        parfor i=1:length(pin)
            pout(i,:) = getMS_local(nirModel, pin(i));
        end
    end
%     It has been handled in build_enir
%     pout = max(pout,0);
%     pout = min(pout,1);
end


function [res] = getMA_local(enirModel, x) 
    SV = enirModel.SV;  % It has already the relative likelihood
    B = length(SV);
    p = zeros(B,1);
    N = length(enirModel.uy);
    % find index of two consecutive training data that x is located among them
    idx = find_idx(enirModel.uy,x); 
    
    assert(length(enirModel.uy)==size(enirModel.beta,1));
    
    beta = enirModel.beta;
    for i=1:B
        if idx == 0
            % We can do extrapolation here
            p(i,:) = beta(1,i);
        elseif idx == N
            % We can do extrapolation here
            p(i,:) = beta(N,i);
        else
            ya = enirModel.uy(idx);
            yb = enirModel.uy(idx+1);
            
            a  = beta(idx,i);
            b  = beta(idx+1,i);
            if (ya == yb)
                p(i,:) = (a+b)/2;
                
                error('Warning: ya==yb: this should not happen. you should handle these cases in the optimization program')

%                 k = idx-1;
%                 tmp = a;
%                 while (k > 0) && (enirModel.y(k) == ya)
%                     tmp = [enirModel.beta(k,i); tmp];
%                     k = k - 1;
%                 end
%                 
%                 k = idx+2;
%                 tmp = [tmp;b];
%                 while (k <= N) && (enirModel.y(k) == yb)
%                     tmp = [tmp;enirModel.beta(k,i)];
%                     k = k + 1;
%                 end
%                 if length(unique(tmp))>1
%                     disp('Ghooghooli ghoo ghoo')
%                 end
%                 p(i,:) = sum(tmp)/length(tmp);  
            else
                p(i,:) = ((x-ya)*b + (yb-x)*a)/(yb-ya);
            end
        end
    end
    res = SV'*p /sum(SV);   
end


function [res] = getMS_local(nirModel, x) 
    N = length(nirModel.y);
    % find index of two consecutive training data that x is located among them
    idx = find_idx(nirModel.y,x); 

    if idx == 0
        % We can do extrapolation here
        res = nirModel.beta(1);
    elseif idx == N
        % We can do extrapolation here
        res = nirModel.beta(N);
    else
        ya = nirModel.y(idx);
        yb = nirModel.y(idx+1);

        a  = nirModel.beta(idx);
        b  = nirModel.beta(idx+1);
        if (ya == yb)
            res = (a+b)/2;  
%             disp('Warning: ya==yb: this should not happen. you should handle these cases in the optimization program')
%             k = idx-1;
%             tmp = a;
%             while (k > 0) && (enirModel.y(k) == ya)
%                 tmp = [enirModel.beta(k,i); tmp];
%                 k = k - 1;
%             end
% 
%             k = idx+2;
%             tmp = [tmp;b];
%             while (k <= N) && (enirModel.y(k) == yb)
%                 tmp = [tmp;enirModel.beta(k,i)];
%                 k = k + 1;
%             end
%             res = sum(tmp)/length(tmp);
        else
            res = ((x-ya)*b + (yb-x)*a)/(yb-ya);
        end
    end   
end




function [res] = find_idx(PTR, x)
    N = length(PTR);
    
    if x>= PTR(end)
        res = N;
        return
    elseif x <= PTR(1)
        res = 0;
        return
    end
    
    minIdx = 1;
    maxIdx = N+1;
    while (maxIdx - minIdx)>1
        midIdx = floor((minIdx+maxIdx)/2);
        if x > PTR(midIdx)
            minIdx = midIdx;
        elseif x < PTR(midIdx)
            maxIdx = midIdx;
        else
            minIdx = midIdx;
            break;
        end
    end

    res = minIdx;
%     Warning: Potentially there could be some odd cased that should be handled
%     here!      
end