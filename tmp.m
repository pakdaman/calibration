
% disp('Using Linear Logistic Regression ')
% [b,dev,model] = glmfit(XTR, YTR, 'binomial', 'link', 'logit');
% [PTR, ~, ~] = glmval(model.beta, XTR, 'logit', model);
% [PTE, ~, ~] = glmval(model.beta, XTE, 'logit', model);








% disp('Using Quadratic Logistic Regression ')
% XTR2 = [XTR, XTR.^2, XTR(:,1).*XTR(:,2), ones(size(XTR,1),1)];
% XTE2 = [XTE, XTE.^2, XTE(:,1).*XTE(:,2), ones(size(XTE,1),1)];
% 
% [b,dev,model] = glmfit(XTR2, YTR, 'binomial', 'link', 'logit');
% [PTR, ~, ~] = glmval(model.beta, XTR2, 'logit', model);
% [PTE, ~, ~] = glmval(model.beta, XTE2, 'logit', model);
% disp('Performance of (Quadratic) LR probabilities is as follows: ')