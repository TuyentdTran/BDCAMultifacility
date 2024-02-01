function W = h2v2( X , A )
    
    [k,d] = size(X); 
    m = size(A,1); 
    
    W = zeros(k,d);
    XAnorm = vecnorm((reshape(X,[],1,d) - reshape(A,1,[],d)),2,3);
    [~,I(:)] = min(XAnorm,[],1);
    for i = 1:k % subtract the\sum_{i=1}^m e_{r(i)} \oprod (x^{r(i)} - a^i) term
        if (any(I==i))
           W(i,:) =  sum((X(i,:)- A(I==i,:))./XAnorm(i,I==i)',1); 
        else
           W(i,:) = 0;
        end
    end
    W = reshape(sum((reshape(X,[],1,d) - reshape(A,1,[],d))./XAnorm,2),[],d) - W; 
end
