function [X,iter,iterlogs]=constrainedBDCA2V2(A,X,projfun,q,parameters,useselfadaptive)
%   A               : (m,n) m data points in R^n.
%   X               : (k,n) initial centers in R^n.
%   projfun         : function handle for projection operator 
%                     that returns a (k,n) array P
%   tau             : Initial penalty parameter
%   sig             : (>1) penalty multiplier.
%   tauf            : (>tau) final penalty.
%   q               : Number of constraints per center
%   useselfadaptive : Flag to turn on self-adaptivity for lambda in BDCA
%   parameters : parameter structure containing the parameters for the
%                algorithm
flag1    = 1; 
iter     = 0;
iterlogs = struct;


tauf = parameters.tauf; 
sig  = parameters.sig; 
mu   = parameters.mu; 
tau  = parameters.tau; 
delta = parameters.delta;
muf   = parameters.muf;

while tau<tauf && flag1 && mu > muf
    iter  = iter+1;
    iterlogs(iter).tau = tau;
    iterlogs(iter).mu  = mu;

    [Xp,bdcaiter,searchiter,lambdaiter,dcanorms] = constrainedBDCAV2(A,X,projfun,tau,mu,q,parameters,useselfadaptive);
    nrm   = norm(Xp-X,'fro')/(norm(X,'fro')+1);
    flag1 = nrm >= parameters.outernorm;
    X     = Xp;
    tau   = sig*tau;
    mu    = delta*mu;
    iterlogs(iter).dcaiter    = bdcaiter; 
    iterlogs(iter).searchiter = searchiter;
    iterlogs(iter).lambdaiter = lambdaiter;
    iterlogs(iter).dcanorms   = dcanorms;
    iterlogs(iter).outernorm  = nrm;
end
end
