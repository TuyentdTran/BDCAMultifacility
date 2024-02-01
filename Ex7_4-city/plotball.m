function plotball(C,r,col,f)
%% plot balls
hold on
kb=size(C,1);
nn=100; t=2*pi/nn*(0:nn); 
e=exp(1i*t); e1=real(e); e2=imag(e);
for j=1:kb
    x=C(j,1)+r(j)*e1; y=C(j,2)+r(j)*e2;
    if nargin<4 || ~f
        plot(x,y,'Color',col,'LineWidth',2)
    else
        h=fill(x,y,col);
        set(h,'EdgeColor','none')
    end
end
hold off
end