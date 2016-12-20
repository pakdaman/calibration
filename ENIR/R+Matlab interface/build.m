function [ enir ] = build( PTR, YTR, bname, option, alpha )
%BUILD_NIR Summary of this function goes here
%   This function builds the near isotonic regression calibration model
%   y: Classifier scores
%   t: the true class of instances (0, 1)
%   bname: the base file name that is used for creating data files
%   communicating R and Matlab

    if nargin < 2
        error('build_nir needs atleast two arguments')
    end
    
    switch nargin
        case 2
            disp('base fie name = data, option = BIC, alpha = 0.005')
            bname = 'data';
            option = 'BIC';
            alpha = 0.005;
        case 3
            disp('option = BIC, alpha = 0.005')
            option = 'BIC';
            alpha = 0.005;
        case 4
            disp('alpha = 0.005')
            alpha = 0.005;
    end
%     Data = [ PTR, YTR];
%     Data = sortrows(Data,[1,2]);
%     y = Data(:,1);
%     z = Data(:,2);
    
    dname = './csv/';
    % Fist delete the existing temporary files
    try
        fname = [dname,bname,'_','z.csv'];
        delete(fname);
    catch
    end
    
    try
        fname = [dname,bname,'_','y.csv'];
        delete(fname);
    catch
    end
    
    try
        fname = [dname,bname,'_','df.csv'];
        delete(fname);
    catch
    end
    
    try
        fname = [dname,bname,'_','lambda.csv'];
        delete(fname);
    catch
    end
    
    try
        fname = [dname,bname,'_','beta.csv'];
        delete(fname);
    catch
    end 
    
    
    
    [y, idx ] = sort(PTR,'ascend');
    z = YTR(idx);
    
    [uy, Iy, Iuy] = unique(y); % Unique y
    %   uy = y(Iy) and y = uy(Iuy).
    
    fname = [dname,bname,'_','z.csv'];
    dlmwrite(fname, z, 'delimiter', ',', 'precision', 50);
    
    fname = [dname,bname,'_','y.csv'];
    dlmwrite(fname, y, 'delimiter', ',', 'precision', 50);
    
    cmd = ['./neariso.R ',' ',dname,' ' ,bname];
    status = -1;
    count = 0;
    maxCount = 10;
    while status ~= 0
        count = count + 1;
        [status,cmdout] = system(cmd);
        if (count > maxCount)
            break;
        end
    end

    if status ==0
        fname = [dname,bname,'_','df.csv'];
        df = csvread(fname);

        fname = [dname,bname,'_','lambda.csv'];
        lambda = csvread(fname);

        fname = [dname,bname,'_','beta.csv'];
        beta = csvread(fname);   
    else
        error('Error in running the R module')
    end
%    Just to Make sure that the predicted values follow rule of probs
    beta = max(beta,0);
    beta = min(beta,1);
    
    % Replacing df with my df which seems to be more reasonable
    df = get_df_local(beta);
    score = zeros(length(lambda),1);
    for i=1:length(lambda)
        p = beta(:,i);
        try
            score(i) = get_score_local(z, p(Iuy), df(i),option);
        catch
            disp('error is here')
        end
    end

    
    idx = ~(isinf(score));
    if length(idx)>1
        idx(1) = (1==0);
    end
%      Remove those that have infinit score as well as the model associated
%      with lambda = 0 (dummy model)
    df = df(idx,:);
    lambda = lambda(idx,:);
    beta = beta(:,idx);
    score = score(idx,:);
    
    MNM = length(df); % Maximum Number of Models
    maxScore = -Inf;
    maxScoreIdx = 0;
    minScore = Inf;
    minScoreIdx = 0;
    SV = zeros(MNM,1);
    for b = 1:MNM
        SV(b,:) = score(b);
        if score(b) > maxScore
            maxScoreIdx = b;
            maxScore = score(b);
        end
        
        if score(b) < minScore
            minScoreIdx = b;
            minScore = score(b);
        end
    end
    SV = exp((min(SV)-SV)./2);% Compute Reletive Likelihood
    
    
    
    model.maxScoreIdx = maxScoreIdx; % Will be used for Model selection
    model.minScoreIdx = minScoreIdx; % Will be used for Model selection
    model.SV = SV;
    model.z = z;
    model.y = y;
    model.uy = uy;
    model.df = df;
    model.beta_org = beta;
    model.lambda = lambda;
    model.Iy = Iy;
    model.Iuy = Iuy;
    
    
    % This part is added to smooth the bin estimates using beta prior like
    % what we did in BBQ 
    model.beta = smooth_beta_local(model);
    
%     WE can do shoulder method as we did for BBQ
    model2 = model;
    if MNM > 1
        idxs = elbow(SV, alpha);
        model2 = processModel_local(model, idxs);
    end    
    
    enir.model = model;
    enir.prunedModel = model2; 
end

function [beta_p] = smooth_beta_local(model)
    % beta_p : process beta (Smoothing probs in the bins using beta prior!)
    z = model.z;
    y = model.y;
%     uy = model.uy;  
    Iy = model.Iy;
    Iuy = model.Iuy; 
    beta = model.beta_org(Iuy,:);
    
    m = size(beta,2); % number of nir models
    n = size(beta,1); % number of training (calibration) instances
%     assert(n == length(uy)); % Assertion to ensure that number of unique y is equal to Beta dimensiton 1
    beta2 = zeros(size(beta));
    for i=1:m
        b = beta(:,i);
        b2 = zeros(size(b));
        j = 1;
        while j <= n
            idx1 = j;
            p_hat = b(idx1);
            j = j + 1;
            while j <= n
                if b(j)~=b(idx1)
                    break;
                end
                j = j + 1;
            end
            idx2 = j - 1;

            o_f = mean(z(idx1:idx2)); % observed frequecny
            e_f = mean(y(idx1:idx2));% expected frequency
            c_e = abs(o_f - e_f); % calibration error in the bin
            cnt = idx2 - idx1 + 1;
            % Smoothing the prob using beta prior(= Laplace Smoothing)
            b2(idx1:idx2) = (p_hat*cnt + e_f *(1-c_e)) / (cnt + 1*(1-c_e));
        end
        beta2(:,i) = b2;
    end
    beta_p = beta2(Iy,:);  
end

function [df] = get_df_local(beta)
    df = zeros(size(beta,2),1);
    for i=1:size(beta,2)
        p = beta(:,i);
        df(i,:) = 1;
        for j=2:length(p)
            if p(j-1)~=p(j)
                df(i,:) = df(i,:) +1;
            end
        end
    end
end
function [ score ] = get_score_local(z, p, df,option)    
    N = length(z);
    l = p;
    idx = (z==0);
    l(idx) = 1-p(idx);
    logLikelihood = sum(log(l));
    
    if strcmp(option,'AICc')
        score = 2*df - 2*logLikelihood + 2*df*(df+1)/(N-df-1);
    elseif strcmp(option,'BIC')
        score = -2*logLikelihood + df *(log(N)+log(2*pi));
%         score = -2*logLikelihood + df *log(N);
    elseif strcmp(option,'AIC')
        score = 2*df -2*logLikelihood;
    else
        error('We only support AIC, AICc, or BIC')
    end    
end


function [outModel] = processModel_local(inModel, idxs)
    outModel.z = inModel.z;
    outModel.y = inModel.y;
    outModel.uy = inModel.uy;
    outModel.SV = inModel.SV(idxs);
    outModel.df = inModel.df(idxs);
    outModel.lambda = inModel.lambda(idxs);
    outModel.beta_org = inModel.beta_org(:,idxs);
    outModel.beta = inModel.beta(:,idxs);
    outModel.maxScoreIdx = 1;
    outModel.minScoreIdx = length(idxs); 
end

%     csvwrite(fname,z);
% function [ score ] = getAICcScore_local( opt )
%     logLikelihood = opt.logLikelihood;
%     N = opt.N;
%     K = opt.K;
%     score = 2*K - 2*logLikelihood + 2*K*(K+1)/(N-K-1);
% end
% 
% function [ score ] = getBICScore_local( opt )
%     logLikelihood = opt.logLikelihood;
%     N = opt.N;
%     K = opt.K;
%     score = -2*logLikelihood + K *(log(N)+log(2*pi));
% end



%             if j == 1
%                 a_i = 0;
%                 b_i = (y(idx2)+y(idx2+1))/2;
%                 m0 = (a_i + b_i) / 2;
%             elseif j == n
%                 a_i = (y(idx1) + y(idx1-1))/2;
%                 b_i = 1;
%                 m0 = (a_i + b_i) / 2;
%             else
%                 a_i = (y(idx1)+y(idx1-1))/2;
%                 b_i = (y(idx2)+y(idx2+1))/2;
%                 m0 = (a_i + b_i) / 2;
%             end

%             p0 = (sum(PTR1)+m0) / (length(PTR1)+1);
%             p_hat2 = (p_hat*cnt + e_f *(1-c_e)) / (cnt + 1*(1-c_e));
%             if p_hat2 > 1 || p_hat2 < 0
%                 error('Havij polo in the build_enir');
%             end