function plotpart(U,Bu,col)
%% plot partion of R^2
hold on
k=size(U,1);
if nargin<3,col=jet(k);
elseif size(col,1)==1,col=repmat(col,k,1);
end
ax=axis;
for j=1:k
    t=[ax(1),ax(2)];
    u=U(j,:);b=Bu(j,:);
    m=-u(1)/u(2);
    if isinf(m)
        y=[ax(3),ax(4)];
    else
        y=m*(t-b(1))+b(2);
        y=[max([y(1),ax(3)]),min([y(2),ax(4)])];
    end
    if m~=0,t=(y-b(2))/m+b(1);end
    if u(1)>0
        x=[t,ax(2),ax(2)];
        if u(2)>0,y=[y,ax(3),ax(4)];
        else,y=[y,ax(4),ax(3)];end
    elseif u(1)<0
        x=[t,ax(1),ax(1)];
        if u(2)>=0,y=[y,ax(4),ax(3)];
        else,y=[y,ax(3),ax(4)];end
    end
    plot(x,y,'Color',col,'LineWidth',1)
end
hold off
end