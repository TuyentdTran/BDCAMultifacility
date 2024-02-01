function P=projpolygon(X,Z)
% returns projection of x onto convex polygon with vertices Z
% X: (k,2) point x in R^2
% Z: {k} cell of k (r,2) matrices representing vertices in R^2
% P: (k,2) nearest point in polygon to x
P=X;


tol = 1e-10;

for j=1:size(X,1)
    x=X(j,:);z=Z{j};
    r = size(z,1);
    % shortcircuit dispatch for rectangles and triangles
    if (r == 3) %triangle 
        % represent x in barycentric coordinates, Ac+z(1,:)' maps
        % from barycentric coordinates to R^2
        A = [z(2,:);z(3,:)]' - z(1,:)'; 
        baryx = A\(x-z(1,:))'; 
        
        %   \|
        %  2 \  3
        %    v3
        %  1 | \   
        %    |  \
        %    |in \
        %-- v1---v2----
        %    |     \
        % 6  |  5   \ 4
        %
        
        logic = [baryx(1)>0,baryx(2)>0,sum(baryx)<1];
        if(isequal([1,1,1],logic))
            %inside do nothing
        else
            %% find nearest vertex
            r=size(z,1);
            y=kron(ones(r,1),x)-z;
            [~,i1]=min(sum(y.^2,2));
            yhat=y(i1,:);
            %% tangent vectors
            t1=z(i1,:)-z(mod(i1-2,r)+1,:);
            t1=t1/norm(t1);
            t2=z(i1,:)-z(mod(i1,r)+1,:);
            t2=t2/norm(t2);
            %% find euclidean projection
            c1=yhat*t1';c2=yhat*t2';
            p=z(i1,:);
            if c1>0 && c2<0
                p=p+c2*t2;
            elseif c1<0 && c2>0
                p=p+c1*t1;
            elseif c1<0 && c2<0
                p1=p+c1*t1;
                p=p+c2*t2;
                [~,l]=min([norm(x-p1),norm(x-p)]);
                if l==1,p=p1;return,end
            end
            P(j,:)=p;
        end
        % rectangle with counter-clockwise ordering of vertices 
    elseif (r == 4 && ( abs(z(1,1)-z(2,1))<tol && abs(z(2,2)-z(3,2))< tol && abs(z(3,1)-z(4,1))<tol && abs(z(4,2)-z(1,2))<tol))
        % handle the 9 possible cases
        
        %  0  | 1  |  2  
        %  ---v3---v2----
        %  3  | in |  4
        % ----v4---v1----
        %  5  |  6 |  7
        
        logic = x<z;
        switch(true)
            case isequal(logical([1,0;1,1;0,1;0,0;]),logic) %in
                % do nothing, keep P(j,:) = X(j,:)
            case isequal(logical([1,0;1,0;1,0;1,0;]),logic) %0
                P(j,:) = z(3,:);
            case isequal(logical([1,0;1,0;0,0;0,0;]),logic) %1
                P(j,:) = [x(1,1), z(2,2)];
            case isequal(logical([0,0;0,0;0,0;0,0;]),logic) %2
                P(j,:) = z(2,:);
            case isequal(logical([1,0;1,1;1,1;1,0;]),logic) %3
                P(j,:) = [z(3,1),x(1,2)];
            case isequal(logical([0,0;0,1;0,1;0,0;]),logic) %4
                P(j,:) = [z(1,1),x(1,2)];
            case isequal(logical([1,1;1,1;1,1;1,1;]),logic) %5
                P(j,:) = z(4,:);
            case isequal(logical([1,1;1,1;0,1;0,1;]),logic) %6
                P(j,:) = [x(1,1), z(1,2)];
            case isequal(logical([0,1;0,1;0,1;0,1;]),logic) %7
                P(j,:) = z(1,:);
            otherwise 
                error("This should not happen");
        end     
    elseif ~inpolygon(x(1),x(2),z(:,1),z(:,2))
        %% find nearest vertex
        r=size(z,1);
        y=kron(ones(r,1),x)-z;
        [~,i1]=min(sum(y.^2,2));
        yhat=y(i1,:);
        %% tangent vectors
        t1=z(i1,:)-z(mod(i1-2,r)+1,:);
        t1=t1/norm(t1);
        t2=z(i1,:)-z(mod(i1,r)+1,:);
        t2=t2/norm(t2);
        %% find euclidean projection
        c1=yhat*t1';c2=yhat*t2';
        p=z(i1,:);
        if c1>0 && c2<0
            p=p+c2*t2;
        elseif c1<0 && c2>0
            p=p+c1*t1;
        elseif c1<0 && c2<0
            p1=p+c1*t1;
            p=p+c2*t2;
            [~,l]=min([norm(x-p1),norm(x-p)]);
            if l==1,p=p1;return,end
        end
        P(j,:)=p;
    end
end
end