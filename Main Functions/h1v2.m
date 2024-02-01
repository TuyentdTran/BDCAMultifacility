function Z = h1v2( X , A ,  mu)
    % evaluation of the Z term in the DCA Multifacility algorithm 
    [k,d] = size(X); 
    m     = size(A,1);
    temp = (reshape(X,[],1,d) - reshape(A,1,[],d))/mu; % all differences of rows between 
    % X - A in the format k x m x d. 

    % subtract unit vectors, part of the ball projection term
    tempnorm = vecnorm(temp,2,3);
    temp = temp - temp./tempnorm;

    % zero out the terms for which temp was already inside the unit ball 
    % and sum 

    temp(reshape(repmat((tempnorm <=1),1,d),k,m,d)) = 0; 
    Z = reshape(sum(temp,2),[],d);
end