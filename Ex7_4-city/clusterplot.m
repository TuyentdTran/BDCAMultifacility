function clusterplot(A,X,I,opt)
%%
if nargin<3,I=[];end
if nargin<4,opt=struct;end
% Default: do not find total center
if ~isfield(opt,'totalCenter'), opt.totalCenter = []; end
% Default data and minor branch options
if ~isfield(opt,'dataStyle'), opt.dataStyle = 'bo'; end
if ~isfield(opt,'dataSize'), opt.dataSize = []; end
if ~isfield(opt,'lineStyle'), opt.lineStyle = 'k'; end
if ~isfield(opt,'lineWidth'), opt.lineWidth = 1; end
% Default centers and major branch options
if ~isfield(opt,'centerStyle'), opt.centerStyle = 'rx'; end
if ~isfield(opt,'centerSize'), opt.centerSize = []; end
if ~isfield(opt,'totalLineStyle'), opt.totalLineStyle = 'k'; end
if ~isfield(opt,'totalLineWidth'), opt.totalLineWidth = 2; end
% Default total center options
if ~isfield(opt,'totalCenterStyle'), opt.totalCenterStyle = 'g*'; end
if ~isfield(opt,'totalSize'), opt.totalSize = []; end
% Default plot options
if ~isfield(opt,'title'), opt.title = []; end
if ~isfield(opt,'pause'), opt.pause = 1/60; end
if ~isfield(opt,'video'), opt.video = []; end

% Number of data points
m=size(A,1);
k=size(X,1);  

% Find indices of nearest centers
if isempty(I), I=center_index(A,X); end
%%
% Plot branches
hold on
for i=1:m
    x=[A(i,1),X(I(i),1)];
    y=[A(i,2),X(I(i),2)];
    plot(x,y,opt.lineStyle,'LineWidth',opt.lineWidth)
end

% Plot total center
if ~isempty(opt.totalCenter)
    x_star = opt.totalCenter;
    for i=1:k
        x=[X(i,1),x_star(1)];
        y=[X(i,2),x_star(2)];
        plot(x,y,opt.totalLineStyle,'LineWidth',opt.totalLineWidth)
    end
    % Plot total center
    b(3)=scatter(x_star(1),x_star(2),opt.totalSize,opt.totalCenterStyle,'filled','MarkerSize',7);
end

% Plot data points and centers
b(1)=scatter(A(:,1),A(:,2),opt.dataSize,opt.dataStyle,'filled');
b(2)=scatter(X(:,1),X(:,2),opt.centerSize,opt.centerStyle,'filled');

% Optional initial centers
if isfield(opt,'X0') 
    if ~isfield(opt,'icenterStyle'), opt.icenterStyle = 'ro'; end
    if ~isfield(opt,'icenterSize'), opt.icenterSize = []; end
    b(length(b)+1)=scatter(opt.X0(:,1),opt.X0(:,2),opt.icenterSize,opt.icenterStyle);    
end
hold off
%%
% Formatting
if isfield(opt,'window'), axis(opt.window); end
if isfield(opt,'legend'), legend(b, opt.legend), end
if isfield(opt,'title'), title(opt.title), end

pause(opt.pause)
if ~isempty(opt.video)
    foo = getframe(gcf);
    writeVideo(opt.video, foo);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function I=center_index(A,X) 
% Find indices of nearest center (if not provided)
    m=size(A,1);
    k=size(X,1);
    I=zeros(1,m);
    for i=1:m
        D=zeros(1,k); % Euclidean distance
        for j = 1:k
            D(j)=norm(A(i,:)-X(j,:));
        end
        % Index of nearest center
        [~,ii]=min(D);
        I(i)=ii(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


