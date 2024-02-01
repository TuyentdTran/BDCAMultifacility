function plotpoly(Z,col)
%% plot ploygons
kp=length(Z);
hold on
for j=1:kp
    z=Z{j};
    plot([z(:,1);z(1,1)],[z(:,2);z(1,2)],'Color',col,'LineWidth',2)
end
hold off
end