function [X,iter,normsiters]=constrainedDCAV2(A,X,projfun,tau,mu,q,parameters)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%              that returns a (k,n) array P
%   tau      : (>0) projection penalty parameter.
%   q        : Number of constraints per center
%   parameters : parameter structure containing the parameters for the
%                algorithm
[k,n] = size(X); 
m     = size(A,1);
S     = ones(k,m)*A;
U     = zeros(k,n);
flag=1;
iter=0;
normsiters = zeros(1000,1); 
while flag
    X_old  = X;
    W      = h2v2(X , A);
    Z      = h1v2(X, A, mu);
    U      = projfun(X);
    Y      = tau*U+W+Z;
    X      = (mu*Y+S)/(m+mu*tau*q);
    
    nrm   = norm(X_old-X,'fro'); 
    normsiters(iter+1) = nrm; 
    iter  = iter+1; 
    
    flag  = nrm >= parameters.innernorm;
end

%remove unneeded entries in the normsiters array
normsiters  = normsiters(1:iter);
end
