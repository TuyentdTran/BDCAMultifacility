function P=projAffineR(X,R,B)
% R: (k,n) jth row orthogonal to hyperplane that partitions R^n
% B: (k,n) jth row is the offset of the hyperplane R+b
% X: (k,n) jth row to be projected onto affine hyperplane

NR = R./vecnorm(R,2,2);
P  = X-(sum((X-B).*NR,2)).*NR; %project onto the hyperplane 
%determine if on the left or right of hyperplane. 
I  = sum((P-X).*NR,2)<0;
P(I,:) = X(I,:);
