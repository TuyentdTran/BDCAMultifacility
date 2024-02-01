function [c,I]=cost(A,X)
%%
m=size(A,1);k=size(X,1);
I=zeros(1,m);c=I;
for i=1:m
    XA=X-repmat(A(i,:),k,1);
    ci=sqrt(sum(XA.^2,2));
    [c(i),I(i)]=min(ci);
end
c=sum(c);
end