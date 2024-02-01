function P=projball(X,C,r)
% returns each projection of X_i onto ball with center C_i and radius r_i
% X: (k,n) k points in R^n
% C: (k,n) centers of k balls
% r: (1,n) radii of k balls
% P: (n,k) P_j is nearest point in Omega_j to X_j
Y = X-C;
N = vecnorm(Y,2,2);  % l-2 norms
% N=norm(Y,2);
r = r./N;
%j=r<1; % indices of vectors inside repective balls
P = r.*Y+C; % project onto the boundary of the ball
inside = r>1;
P(inside,:) = X(inside,:);  % if X is inside the inside(i) ball, keep the same value
end