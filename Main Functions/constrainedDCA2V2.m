function [X,iter,iterlogs]=constrainedDCA2V2(A,X,projfun,q,parameters)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%              that returns a (k,n) array P
%   q        : Number of constraints per center
%   parameters : parameter structure containing the parameters for the
%                algorithm
flag1    = 1;
iterlogs = struct;
iter     = 0;

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

    [Xp,dcaiter,dcanorms] = constrainedDCAV2(A,X,projfun,tau,mu,q,parameters);
    nrm   = norm(Xp-X,'fro')/(norm(X,'fro')+1);
    flag1 = nrm >= parameters.outernorm;
    X     = Xp;
    tau   = sig*tau;
    mu    = delta*mu;

    iterlogs(iter).dcaiter    = dcaiter;
    iterlogs(iter).outernorms = nrm;
    iterlogs(iter).dcanorms   = dcanorms;
end