%%%%
% new_76 - Test DCA,  BDCA adaptive for a MLP on the EIL76 data set. This is used to generate 
% a single run for each, and then plot the resulting solution.
% 
% See example 2 in the paper. 
%
%%%%


% Setup a parameter structure for the algorithm

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
parameters.lambdahistory = 1;
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
projfun=@(mat)proj(mat,C,r,Z,q);

%initial center
a=2*pi*rand(1,2);
r_rand=sqrt(rand(1,2));

sel = [1,4];
X0  = (r(sel).*r_rand).*[ cos(a); sin(a)]+C(sel,:); 

C2 = [32.8 37; 40 18.2];
r2 = [2.7,1.65];

tic
[X,N1,iterlogsdca]        = constrainedDCA2V2(A,X0,projfun,q,parameters);
toc
[c,I]=cost(A,X);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA= %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%% BDCA Adaptive
% profile on
tic
[X3,outeriter,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,q,parameters,true);
toc
[c3,I3]=cost(A,X3);
bdcaiter_adap = 0;
for i = 1:outeriter
    bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA adaptive = %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter_adap,c3);


%% plots (we can use X or X3 here)
fig=figure();clf
hold on
wid=2;
axis equal

for j=1:5
    draw_circle(C(j,1),C(j,2),r(j),'k',[]);
end
plotpoly(Z,'k')

plotclusters(A,X3,I3)
% plot(X0(:,1),X0(:,2),'o');

hold off
grid on

%% Projection function
function P=proj(X,C,r,Z,q)
[~,n]=size(X);
Xt = [repmat(X(1,:),2,1); repmat(X(2,:),3,1);];
Pb = projball(Xt,C,r);
Pp=projpolygon(X(1,:),Z);
    T  = [Pp; Pb];
    P  = reshape(sum(reshape(T, q, [])), [], size(T, 2));   
end