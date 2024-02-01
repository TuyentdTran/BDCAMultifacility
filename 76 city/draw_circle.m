% Draw a circle in R2
function draw_circle (xcenter,ycenter,radius,color,option)
format long
theta = linspace(0,2*pi,200);
x = xcenter + radius * cos(theta);
y = ycenter + radius * sin(theta);
    if option==true
        plot(x,y,'Color',color,'Linewidth',0.7,'Linestyle','--');
    else
        plot(x,y,'Color',color,'Linewidth',1.5);
    end
end
