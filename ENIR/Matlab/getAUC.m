function [AUC] = getAUC (Actual,Predicted)
if length(unique(Actual))~=2 ||max(unique(Actual))~=1
    error('Strange input for AUC analysis');
end
nTarget     = sum(double(Actual == 1));
nBackground = sum(double(Actual ~= 1));

% Rank data
R = tiedrank(Predicted);  % 'tiedrank' from Statistics Toolbox

% Calculate AUC
AUC = (sum(R(Actual == 1)) - (nTarget^2 + nTarget)/2) / (nTarget * nBackground);

AUC = max(AUC,1-AUC);