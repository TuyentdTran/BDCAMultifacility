function [LogsDCA,LogsBDCA_adapt,parameters] = loopmain(N)
% Run a test N times
% parameters for the solver 
% Setup a parameter structure for the algorithm

parameters = struct; 
parameters.tau = 1;  % projection penalty
parameters.sig = 100;% penalty growth parameter
parameters.tauf = 1e8; 
parameters.mu = 1; %smoothing parameter
parameters.delta = 0.25;
parameters.muf   = 1e-6;

parameters.outernorm = 1e-6;
parameters.innernorm = 1e-6;    

% BDCA specific parameters

parameters.alp  = 0.05;
parameters.beta = 0.01;
parameters.lambdastart = 1;
parameters.gamma = 2;
parameters.lambdahistory = 1;
parameters.maxsearch = 2; 
parameters.lambdamin = 1e-3;
parameters.lambda_skip = 100;  


%% load data set
A = dlmread('continental.dat');
[m,n]=size(A);

%% projection function

k=4; % number of centers
q=2;

% centers within r(i) degrees lat/long of cities, i=1,..,4
C=[	-116.6873596	43.6629384;%Caldwell,Idaho
    -97.54021069 35.4715076;% Oklahoma City, Oklahoma
    -87.7333934  42.0333636;% Skokie,Illinois
     -73.98143112	40.76149276; %New York,New York
     -77.0363658	38.8951118]; %Washington,District of Columbia
r=[4 6 2 4 4]';

% centers lying in polygons
Z={[-102.05,37;-102.05,41;-109.05,41;-109.05,37],...% state of colorado
    [-93.6091064,41.6005448;...% des moine
    -81.6954088,41.4994954;...% cleveland
    -84.3879824,33.7489954]};% atlanta
% kp=length(Z);
% center lying west of -115 lon
U  = [-1,0];
Bu = [-115,0];
% kpart=size(U,1);

projfun = @(mat)projcities(mat,C,r,Z,U,Bu,q);

for runiter = 1:N          
    % Starting point
    X0= zeros(k,n);
    for i = 1:k-1
        X0(i,:)=C(i,:)+r(i)*rand(1,2)/sqrt(2);
    end
    X0(4,:)=C(5,:)+r(5)*rand(1,2)/sqrt(2);
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
function P=projcities(X,C,r,Z,U,Bu,q)
Y=[X;X(4,:)];
Pb=projball(Y,C,r);
Q=X(2:3,:);
Pp=projpolygon(Q,Z);
Pt=projAffineR(X(1,:),U,Bu);
T = [Pt;Pb(1,:);Pb(2,:);Pp(1,:);Pb(3,:);Pp(2,:);Pb(4,:);Pb(5,:)];
P  = reshape(sum(reshape(T, q, [])), [], size(T, 2));    
end