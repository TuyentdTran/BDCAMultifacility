function P=projball2(X,C,r)
    % returns each projection of a X_i onto ball with center C_i and radius r_i
    % X: (1,n) point in R^n
    % C: (k,n) centers of k balls
    % r: (1,n) radii of k balls
    % P: (k,n) P_j is nearest point in Omega_j to X_j
    if(size(X,1) > 1) 
        error("X must be single point in R^n"); 
    end
    Y = X-C;
    N = vecnorm(Y,2,2); %sum(Y.^2,2).^(1/2); % l-2 norms
    r = r./N';
    
    P = r.*Y+C; % project onto the boundary of the ball
    inside = r>1;
    P(inside,:) = X(inside,:);  % if X is inside the inside(i) ball, keep the same value
end