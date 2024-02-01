function [X,iter,searchiters,lambdaiters,normsiters]=constrainedBDCAV2_Test(A,X,projfun,tau,mu,q,parameters,useselfadaptive)
%   A               : (m,n) m data points in R^n.
%   X               : (k,n) initial centers in R^n.
%   projfun         : function handle for projection operator 
%                     that returns a (k,n) array P
%   tau             : (>0) projection penalty parameter.
%   q               : Number of constraints per center
%   useselfadaptive : Flag to turn on self-adaptivity for lambda in BDCA 
%   parameters : parameter structure containing the parameters for the
%                algorithm

%   alpha,beta and lambda are parameters that will be used in BDCA and set
%   inside of the function.

[k,n]       = size(X); 
m           = size(A,1);
alp         = parameters.alp; 
beta        = parameters.beta;
lambdastart = parameters.lambdastart; 
gamma       = parameters.gamma;
lambdahistorylength = parameters.lambdahistory;
max_search = parameters.maxsearch;
lambda_min = parameters.lambdamin;
lambda_skip = parameters.lambda_skip;

S=ones(k,m)*A;
flag1 = 1;
iter  = 0;
searchiters  = zeros(1000,1);
lambdaiters  = zeros(1000,1);
normsiters   = zeros(1000,1);
lambdahist   = zeros(lambdahistorylength,1);
lambdabar    = zeros(lambdahistorylength,1);
lambda       = lambdastart;
lambda_skip_iter = lambda_skip; 

while flag1
    X_old  = X;
    W      = h2v2(X , A);
    Z      = h1v2(X, A, mu);
    U      = projfun(X);
    Y      = tau*U+W+Z;
    Xk     = (mu*Y+S)/(m+mu*tau*q);
    
    if lambda < lambda_min && lambda_skip_iter == lambda_skip
        % skip line search for the next lambda_skip steps. 
        lambdahist = lambda_min; 
        lambda = gamma*lambda_min;
        searchiters(iter+1) = -1;
        lambdaiters(iter+1) = 0;
        lambda_skip_iter = 0; 
        X = Xk;
    elseif lambda_skip_iter < lambda_skip
        % skip line search
        lambda_skip_iter = lambda_skip_iter + 1; 
        searchiters(iter+1) = -1;
        lambdaiters(iter+1) = 0;
        X = Xk;
    else
        if(useselfadaptive)
            %renewhistory
            lambdahist(1:end-1) = lambdahist(2:end);
            lambdahist(end)     = lambda;

            if(all(abs(lambdabar - lambdahist)<1e-6))
                %no search was done for the last lambdahistorylength iterations
                %increase lambda from previous        
                lambda = gamma*lambdahist(end); 
            else % use previous lambda
                lambda = lambdahist(end);
            end
            %renew history 
            lambdabar(1:end-1)  = lambdabar(2:end);
            lambdabar(end)      = lambda;
        else
            lambda = lambdastart;
        end
        d_k        = Xk - X_old;
        X          = Xk + (d_k * lambda);
        phi_k      = cost_pen_opt(A,Xk,tau,q,mu,projfun);
        costpen    = cost_pen_opt(A,X,tau,q,mu,projfun); 
        rhs        = alp * norm(d_k, 'fro')^2;
        d          = phi_k - lambda^2 * rhs;
        searchiter = 0; 
        while costpen > d && searchiter < max_search
                lambda     = lambda * beta;
                X          = Xk + (d_k * lambda);
                costpen    = cost_pen_opt(A,X,tau,q,mu,projfun);
                d          = phi_k - lambda^2 * rhs;
                searchiter = searchiter + 1;
        end
        searchiters(iter+1) = searchiter;
        lambdaiters(iter+1) = lambda;
    end
      
    nrm  = norm(X_old-X,'fro');
    normsiters(iter+1) = nrm;
    iter = iter+1;
    
    flag1 = nrm >= parameters.innernorm;
end
%remove unneeded entries in the searchiters array
searchiters = searchiters(1:iter);
lambdaiters = lambdaiters(1:iter);
normsiters  = normsiters(1:iter);
end