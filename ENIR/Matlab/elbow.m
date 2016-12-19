function [res] = elbow(scores, alpha)
    b = length(scores);
    sigma2 = var(scores,1);
    [R, idxs ] = sort(scores,'descend');
    k = 1;
    while (R(k)==R(k+1))
        k = k + 1;
    end
    while (k<b)&&((R(k)-R(k+1))/sigma2 > alpha)
        k = k+1;
    end
    
    if k  > 1
        res = idxs(1:k-1);
    else
        res = idxs(1);
    end
    
end