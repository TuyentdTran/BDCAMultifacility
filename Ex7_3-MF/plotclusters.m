function plotclusters(A,X,I)
% Number of data points
m=size(A,1); 

%%
% Plot to connect data points and centers
hold on
for i=1:m
    x=[A(i,1),X(I(i),1)];
    y=[A(i,2),X(I(i),2)];
    plot(x,y,'magenta','Linewidth',0.7)
end
%Plot data points
scatter(A(:,1),A(:,2),10,'bo','filled');
% Plot centers
plot(X(:,1),X(:,2),'bp','MarkerFaceColor','b','MarkerSize',20)
end