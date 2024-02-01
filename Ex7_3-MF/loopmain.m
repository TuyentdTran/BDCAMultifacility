function [LogsDCA,LogsBDCA_adapt,parameters] = loopmain(N)
% Run a test N times

% parameters for the solver 
% Setup a parameter structure for the algorithm

parameters = struct; 
parameters.tau = 1;  % penalty parameter
parameters.sig = 10;% penalty growth parameter
parameters.tauf = 1e8; 
parameters.mu = 1; %smoothing parameter
parameters.delta = 0.25; %0.75;
parameters.muf = 1e-6;

parameters.outernorm = 1e-6;
parameters.innernorm = 1e-6; 

% BDCA specific parameters

parameters.alp  = 0.05;
parameters.beta = 0.01;
parameters.lambdastart = 0.1;
parameters.gamma = 2;
parameters.lambdahistory = 2;
parameters.maxsearch = 2; 
parameters.lambdamin = 1e-3;
parameters.lambda_skip = 60; 

%% load data set

x1=[2 4 4 2]; 
y1=[2 2 4 4];
rc=0.3*ones(4,1);

m=1000;
A=zeros(m,2);


for j=1:4
    for t=1:(m/4) %loop until doing n points inside the circle
        if j>1
            s=(m/4)*(j-1)+t;
        else
            s=t;
        end
    [x,y]=cirrdnPJ(x1(j),y1(j),rc(j));
    A(s,:)=[x,y];
    end
end
[m,n]=size(A);
%% projection function
C=[3 3;3 3;3 3;3 3]; % ball centers
q=1; 
r=0.3*ones(4,1);

projfun=@(mat)projball(mat,C,r);


for runiter = 1:N          
%     X0 = rand(4,2)/sqrt(2).*rc + [x1',y1'];
    X0 = rand(4,2).*rc + [x1',y1'];
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

% %% Projection function
% function P=projhardcoded(X,C,r,Z)
%     [~,n]=size(X);
%     Pb=zeros(3,n);
%     Pb(1,:)=projball2(X(1,:),C(1,:),r(1));
%     Pb(2,:)=projball2(X(2,:),C(2,:),r(2));
%     Pb(3,:)=projball2(X(2,:),C(3,:),r(3));
%     Pp=projpolygon(X(1,:),Z);
%     P=[Pp(1,:);Pb(1,:);Pb(2,:);Pb(3,:)];
% end

