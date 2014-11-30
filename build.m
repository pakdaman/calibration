function [ bbq ] = build( PTR, YTR, options )
% Input:
% PTR: Vector of predicted values
% YTR: Vector of true labels {0,1}
% Options: Structure of optional choices
       
    if isfield(options,'scoringFunc')
        scoringFunc = options.scoringFunc;
    else
        scoringFunc = 'BDeu2';
        % Special case of BDeu inwhich N0 = 2B where B is number of bins
    end
    
    
    % Used in case that the scoring function is BDeu (we need to set N0)
    if isfield(options,'N0')
        N0 = options.N0;
    else
        N0 = 2;
    end
    
    
    % The threshold parameter for pruning final models based on their score
    % using elbow method
    if isfield(options,'alpha')
        alpha = options.alpha;
    else
        alpha = 0.001;
    end
    
    
    if isfield(options,'sort') 
        % Useful for prevent sorting if the predictions are already sorted
        runSort = options.sort;
    else
        runSort = 1;
    end
    
    
    if runSort == 1
        [~, idx] = sort(PTR);
        PTR = PTR(idx);
        YTR = YTR(idx);
    end
        
    N = length(PTR);
    lnfact = initlnfact_local(N + 1); % array includes log(n!) : starting from 0 ... N 
    % i.e.  lnfact(1) = 0!, lnfact(2) = 1!, lnfact(3) = 2!; ...
    
    maxBinNo = min(ceil(N/5),ceil(10*N^(1/3))); % To have at least 5 instances in each bin !
    minBinNo = max(1, floor(N^(1/3)/10)); % To have at lease 1 bin !
    MNM = maxBinNo-minBinNo +1; % maximum number of models
    model = cell(MNM,1);
    model{1}.scoringFunc = scoringFunc;
    
    opt1.PTR = PTR;
    opt1.lnfact = lnfact;
    opt1.N0 = N0;
    if strcmp(scoringFunc,'K2')
        scoreFunc = @getK2Score_local;
    elseif strcmp(scoringFunc,'AIC')
        scoreFunc = @getAICScore_local;
    elseif strcmp(scoringFunc,'AICc')
        scoreFunc = @getAICcScore_local;
    elseif strcmp(scoringFunc,'BIC')
        scoreFunc = @getBICScore_local;
    elseif strcmp(scoringFunc,'BDeu') % N0 = 2 by default
        scoreFunc = @getBDeuScore_local;
    elseif strcmp(scoringFunc,'BDeu2') % N0 = 2B by default (the version that used in the paper)
        scoreFunc = @getBDeuScore2_local;    
    end
    
    
    parfor b = 1 : MNM 
        [histModel, cutIdx, cutPoints, logLikelihood]  = HistCalibrationFreq_local( PTR, YTR , b + minBinNo -1);
        funcOpt = buildFuncOpt_local(opt1, histModel, cutIdx, cutPoints, logLikelihood);
        score = scoreFunc(funcOpt);
        model{b}.histModel= histModel;
        model{b}.cutIdx = cutIdx;
        model{b}.cutPoints = cutPoints;
        model{b}.score = score;
        model{b}.logLikelihood = logLikelihood;
    end
    
    maxScore = -Inf;
    maxScoreIdx = 0;
    minScore = Inf;
    minScoreIdx = 0;
    SV = zeros(MNM,1);
    for b = 1:MNM
        SV(b,:) = model{b}.score;
        if model{b}.score > maxScore
            maxScoreIdx = b;
            maxScore = model{b}.score;
        end
        
        if model{b}.score < minScore
            minScoreIdx = b;
            minScore = model{b}.score;
        end
    end
    SV = exp((min(SV)-SV)./2);% Compute Reletive Likelihood
    
    model{1}.maxScoreIdx = maxScoreIdx; % Will be used for Model selection
    model{1}.minScoreIdx = minScoreIdx; % Will be used for Model selection
    model{1}.SV = SV;
    idxs = elbow(SV, alpha);
    model2 = processModel_local(model, idxs);
    bbq.model = model;
    bbq.prunedModel = model2;
end

function [outModel] = processModel_local(inModel, idxs)
    for i = 1:length(idxs)
        outModel{i} = inModel{idxs(i)};
    end
    outModel{1}.minScoreIdx = 1;
    outModel{1}.SV = inModel{1}.SV(idxs);
end


function [funcOpt] = buildFuncOpt_local(opt, histModel, cutIdx, cutPoints, logLikelihood)
    N = length(opt.PTR);
    K = length(histModel);
    
    funcOpt.histModel = histModel;
    funcOpt.cutIdx = cutIdx;
    funcOpt.cutPoints = cutPoints;
    funcOpt.logLikelihood = logLikelihood;
    funcOpt.K = K;
    funcOpt.N = N;
    funcOpt.PTR = opt.PTR;
    funcOpt.lnfact = opt.lnfact;
    funcOpt.N0 = opt.N0;
end

function [ score ] = getBDeuScore_local( opt )   
    histModel = opt.histModel;
 
    % Firs Compute lnMarginal
    B = length(histModel);
    N0 = opt.N0;
    C = N0/B;
    
    score = B*gammaln(C);
    for j=1:B
        nj = histModel{j}.n;
        nj0 = histModel{j}.n0;
        nj1 = histModel{j}.n1;
        pj = 0.5;
        score = score + gammaln(nj1+C*pj)+gammaln(nj0+C*(1-pj)) - gammaln(nj+C) ...
                      - gammaln(C*pj) - gammaln(C*(1-pj)); 
    end
    
% It is possible to add the effect of structural prior here, however we used uniform
% prior for now
%     x = log(B/(opt.N^(1/3)));
%     score = score + log(normpdf(x,0,0.75));
% end of adding prior
    score = -2*score;
end

function [ score ] = getAICcScore_local( opt )
    logLikelihood = opt.logLikelihood;
    N = opt.N;
    K = opt.K;
    score = 2*K - 2*logLikelihood + 2*K*(K+1)/(N-K-1);
end


function [ score ] = getBICScore_local( opt )
    logLikelihood = opt.logLikelihood;
    N = opt.N;
    K = opt.K;
    score = -2*logLikelihood + K *(log(N)+log(2*pi));
end


function [ score ] = getK2Score_local( opt )
    histModel = opt.histModel;
    lnfact = opt.lnfact;
    
    score = 0;
    % Firs Compute lnMarginal
    B = length(histModel);
    
    for b=1:B
        U = [histModel{b}.n0, histModel{b}.n1];
        score = score + lnMarginalLikelihood_local(U , lnfact);
    end
    
% I can add the effect of structural prior here, however I use uniform
% prior for now
%     x = log(B/(opt.N^(1/3)));
%     score = score + log(normpdf(x,0,0.75));
% end of adding prior

    score = -2*score;
end


function [ score ] = getBDeuScore2_local( opt )
    histModel = opt.histModel;
   
    % Firs Compute lnMarginal
    B = length(histModel);
    N0 = 2*B;
    C = N0/B;
    
    score = B*gammaln(C);
    for j=1:B
        nj = histModel{j}.n;
        nj0 = histModel{j}.n0;
        nj1 = histModel{j}.n1;
        pj = (histModel{j}.min + histModel{j}.max)/2;
        pj = min(pj,1 - 5*10^-3);
        pj = max(pj, 5*10^-3);
        score = score + gammaln(nj1+C*pj)+gammaln(nj0+C*(1-pj)) - gammaln(nj+C) ...
                      - gammaln(C*pj) - gammaln(C*(1-pj)); 
    end
    
% I can add the effect of structural prior here, however I use uniform
% prior for now
%     x = log(B/(opt.N^(1/3)));
%     score = score + log(normpdf(x,0,0.75));
% end of adding prior
    score = -2*score;
end


function [res] = lnMarginalLikelihood_local(U, lnfact)
    tmp = 0;
    pad = 1; % Because Matlab Array start from 1 not 0
    for i=1:length(U)
        tmp = tmp + lnfact(U(i)+pad);
    end
    res = lnfact(1+pad) - lnfact(1+sum(U)+pad) + tmp;
end



function [lnfact] = initlnfact_local(N)
  lnfact = zeros(N+1,1);
  for w = 2:N+1
      lnfact(w) = lnfact(w-1) + log(w-1);
  end
end


function [ histModel, cutIdx, cutPoints, logLikelihood ] = HistCalibrationFreq_local( PTR, YTR , b)
    %HISTCALIBRATION Summary of this function goes here
    %   Detailed explanation goes here

    % Input Invariant: I assume that PTR is sorted
    % compute logMarginalLikelihood based on BDeu scoring measure
    
    N = length(YTR);
    
    logLikelihood = 0;
    if b==1
        histModel{b}.min = 0;
        histModel{b}.max = 1; 
        % Building an intuitive prior for Beta smoothing
        m0 = (histModel{b}.min + histModel{b}.max)/2;
        idx = (YTR==1);
        PTR1 = PTR(idx);
        p0 = (sum(PTR1)+m0) / (length(PTR1)+1);
        % end of building intuitive prior
        
        histModel{b}.n = length(YTR);
        histModel{b}.n1 = sum(YTR);
        histModel{b}.n0 = histModel{b}.n - histModel{b}.n1;
        histModel{b}.P = (histModel{b}.n1+p0) / (histModel{b}.n+1); 
        
        if histModel{b}.n1 > 0
            logLikelihood = logLikelihood + histModel{b}.n1*log(histModel{b}.P);
        end
        if histModel{b}.n0 > 0
            logLikelihood = logLikelihood + histModel{b}.n0*log(1-histModel{b}.P);
        end
        
        cutIdx = [];
        cutPoints = [];
        
        
    else % Where b > 1
               
        maxNj = 0;
        Yhat = PTR;
        Y = YTR;
        c = floor(length(Y)/b);
        i = 1;
        idx = 1;

        while (i < b)
            idx1 = (i-1)*c + 1;
            idx2 = i*c;
            j = i + 1;
            while j <= b
                if j<b
                    Jidx2 = j*c;
                    if PTR(Jidx2) == PTR(idx1)
                        idx2 = Jidx2;
                        j = j + 1;    
                    else
                        break;
                    end
                else
                    Jidx2 = N;
                    if PTR(Jidx2) == PTR(idx1)
                        idx2 = Jidx2;
                        j = j + 1;    
                    else
                        break;
                    end
                end
            end
            
            T{idx}.Y = Y(idx1:idx2,:);
            T{idx}.PTR = PTR(idx1:idx2,:);
            T{idx}.Yhat = Yhat(idx1:idx2,:);
            maxNj = max(maxNj , idx2-idx1+1);
            if idx2 < N
                cutIdx(idx) = idx2; 
            end
            idx = idx + 1;
            i = j;
        end
        
        if idx2 < N
            T{idx}.Y = Y(idx2+1:end,:);
            T{idx}.PTR = PTR(idx2+1:end,:);
            T{idx}.Yhat = Yhat(idx2+1:end,:);
        end
        
        
        b0 = b;
        b = length(T);  
        histModel = cell(length(T),1);

        histModel{1}.min = 0;
        histModel{1}.max = (T{1}.Yhat(end)+T{2}.Yhat(1))/2;
        cutPoints(1,:) = histModel{1}.max;
        
        
        % Building an intuitive prior for Beta smoothing
        m0 = (histModel{1}.min + histModel{1}.max)/2;
        idx = (T{1}.Y==1);
        PTR1 = T{1}.PTR(idx);
        p0 = (sum(PTR1)+m0) / (length(PTR1)+1);
        % Enf of building prior
        
        histModel{1}.n = length(T{1}.Y);
        histModel{1}.n1 = sum(T{1}.Y);
        histModel{1}.n0 = histModel{1}.n - histModel{1}.n1;
        histModel{1}.P = (histModel{1}.n1+p0) / (histModel{1}.n+1); 
        if histModel{1}.n1 > 0
            logLikelihood = logLikelihood + histModel{1}.n1*log(histModel{1}.P);
        end
        if histModel{1}.n0 > 0
            logLikelihood = logLikelihood + histModel{1}.n0*log(1-histModel{1}.P);
        end
        
                
        for i=2:b-1
            histModel{i}.min = (T{i}.Yhat(1)+T{i-1}.Yhat(end))/2;
            histModel{i}.max = (T{i}.Yhat(end)+T{i+1}.Yhat(1))/2;
            cutPoints(i,:) = histModel{i}.max;
            
            
            % Building an intuitive prior for Beta smoothing
            m0 = (histModel{i}.min + histModel{i}.max)/2;
            idx = (T{i}.Y==1);
            PTR1 = T{i}.PTR(idx);
            p0 = (sum(PTR1)+m0) / (length(PTR1)+1);
            % Enf of building prior
            
            histModel{i}.n = length(T{i}.Y);
            histModel{i}.n1 = sum(T{i}.Y);
            histModel{i}.n0 = histModel{i}.n - histModel{i}.n1;
            histModel{i}.P = (histModel{i}.n1+p0) / (histModel{i}.n+1); 
            if histModel{i}.n1 > 0
                logLikelihood = logLikelihood + histModel{i}.n1*log(histModel{i}.P);
            end
            if histModel{i}.n0 > 0
                logLikelihood = logLikelihood + histModel{i}.n0*log(1-histModel{i}.P);
            end
            
        end
        histModel{b}.min = (T{b}.Yhat(1)+T{b-1}.Yhat(end))/2;
        histModel{b}.max = 1; 
    
        % Building an intuitive prior for Beta smoothing
        m0 = (histModel{b}.min + histModel{b}.max)/2;
        idx = (T{b}.Y==1);
        PTR1 = T{b}.PTR(idx);
        p0 = (sum(PTR1)+m0) / (length(PTR1)+1);
        % Enf of building prior
        
        
        histModel{b}.n = length(T{b}.Y);
        histModel{b}.n1 = sum(T{b}.Y);
        histModel{b}.n0 = histModel{b}.n - histModel{b}.n1;
        histModel{b}.P = (histModel{b}.n1+p0) / (histModel{b}.n + 1); 
        if histModel{b}.n1 > 0
            logLikelihood = logLikelihood + histModel{b}.n1*log(histModel{b}.P);
        end
        if histModel{b}.n0 > 0
            logLikelihood = logLikelihood + histModel{b}.n0*log(1-histModel{b}.P);
        end
        
    end

end

