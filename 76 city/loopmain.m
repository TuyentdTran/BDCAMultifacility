function [LogsDCA,LogsBDCA_adapt,parameters] = loopmain(N)
% Run a test N times

parameters = struct; 
parameters.tau = 1;  % projection penalty
parameters.sig = 1e2;% penalty growth parameter
parameters.tauf = 1e8; 
parameters.mu = 1; %smoothing parameter
parameters.delta = 0.85;
parameters.muf   = 1e-6;

parameters.outernorm = 1e-6;
parameters.innernorm = 1e-6;    

% BDCA specific parameters

parameters.alp  = 0.05;
parameters.beta = 1e-2;
parameters.lambdastart = 1;
parameters.gamma = 2;
parameters.lambdahistory = 2;
parameters.maxsearch = 2; 
parameters.lambdamin = 1e-3;
parameters.lambda_skip = 150;  

%% load data set
A = dlmread('eil76.dat');
%% projection function
k = 2; %number of constraints per center, assumed 2 elsewhere
q = 3;
%constraint balls
C     = [30 40; 32 37; 35 20; 45 20; 40 15];

r = [7,4.5,7,7,5]';
% center lying in polygon
Z  = {[40,30; 40,40; 30,40; 30,30]};

projfun = @(mat) proj(mat,C,r,Z,q);

for runiter = 1:N          
    X0=[C(1,:)+r(1)*rand(1,2)/sqrt(2) ;C(4,:)+r(4)*rand(1,2)/sqrt(2)]; 
    %% DCA
    tic
    [X,N1,iterlogsdca]        = constrainedDCA2V2(A,X0,projfun,q,parameters);
    LogsDCA(N-runiter+1).time = toc;
    LogsDCA(N-runiter+1).logs = iterlogsdca;
    c = cost(A,X);
    dcaiter = 0;
    for i = 1:N1
        dcaiter = dcaiter + iterlogsdca(i).dcaiter;
    end
    LogsDCA(N-runiter+1).totaliter = dcaiter;
    LogsDCA(N-runiter+1).outeriter = N1;
    LogsDCA(N-runiter+1).cost      = c;
    %% BDCA -- Adaptive
    tic
    [X2,outeriter,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,q,parameters,true);
    LogsBDCA_adapt(N-runiter+1).time = toc;
    LogsBDCA_adapt(N-runiter+1).logs      = iterlogsbdca_adap;
    c2 = cost(A,X2);
    bdcaiter_adap = 0;
    for i = 1:outeriter
        bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
    end
    LogsBDCA_adapt(N-runiter+1).totaliter = bdcaiter_adap; 
    LogsBDCA_adapt(N-runiter+1).outeriter = outeriter; 
    LogsBDCA_adapt(N-runiter+1).cost      = c2;
end

end

%% Projection function
function P=proj(X,C,r,Z,q)
[~,n]=size(X);
Xt = [repmat(X(1,:),2,1); repmat(X(2,:),3,1);];
Pb = projball(Xt,C,r);
Pp=projpolygon(X(1,:),Z);
    T  = [Pp; Pb];
    P  = reshape(sum(reshape(T, q, [])), [], size(T, 2));   
end
