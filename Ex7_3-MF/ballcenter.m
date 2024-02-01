clear;
format long

% parameters for the solver 
% Setup a parameter structure for the algorithm

parameters = struct; 
parameters.tau = 1;  % projection penalty
parameters.sig = 10;% penalty growth parameter
parameters.tauf = 1e8; 
parameters.mu = 1; %smoothing parameter
parameters.delta = 0.25;
parameters.muf = 1e-6;

parameters.outernorm = 1e-6;
parameters.innernorm = 1e-6; 

% BDCA specific parameters

parameters.alp  = 0.05;
parameters.beta = 0.01;
parameters.lambdastart = 0.1;
parameters.gamma = 2;
parameters.lambdahistory = 1;
parameters.maxsearch = 2; 
parameters.lambdamin = 1e-3;
parameters.lambda_skip = 70;  

% Problem Setup 
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
   % plot(x,y,'bx');
   % pause(0.1)
    end
end

[m,n]=size(A);

%% projection function
C=[3 3;3 3;3 3;3 3]; % ball centers
q=1; 
k = 4;%number of centers
r=0.3*ones(k,1);


projfun=@(mat)projball(mat,C,r);

%% initializaions
X0 = rand(k,2).*rc + [x1;y1]';

%% DCA
% profile on
tic
[X,N1,iterlogsdca]=constrainedDCA2V2(A,X0,projfun,q,parameters);
toc
% dca_prof = profile("info");
[c,I]=cost(A,X);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA= %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%% BDCA Adaptive
% profile on
tic
[X3,outeriter_adap,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,q,parameters,true);
toc
% bdca_adapt_prof=profile("info");
[c3,I3]=cost(A,X3);
bdcaiter_adap = 0;
for i = 1:outeriter_adap
    bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA adaptive = %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter_adap,bdcaiter_adap,c3);

% plots (we can use X or X3 here)
fig=figure();
axis equal
hold on


for j=1:4
    draw_circle(x1(j),y1(j),rc(j),'k',[]);
end
hold on
draw_circle(C(1,1),C(1,2),r(1),'b',[])
hold on
plotclusters(A,X3,I3)
hold off
% grid on