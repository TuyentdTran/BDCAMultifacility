function c=cost_pen_opt(A,X,tau,q,mu,projfun)

% An optimized cost penalty function that doesn't use DCA splitting to 
% evaluate 

d=size(A,2);

temp = (reshape(X,[],1,d) - reshape(A,1,[],d))/mu; % all differences of rows between 
% X - A in the format k x m x d. 

% subtract unit vectors, part of the ball projection term
XAnorm  = vecnorm(temp,2,3);
disnorm = vecnorm(temp - temp./XAnorm,2,3).^2;

disnorm(XAnorm<=1) = 0;
 
T = projfun(X);
c = mu/2*sum(min(XAnorm.^2-disnorm,[],2))+tau/2*(sum(vecnorm(X-T,2,2).^2)); % this 2nd term is hardcoded
end