clear all
clc
clf

%% forest
x = 0:1:1;
y = 0:1:1;

r(length(x)-1, length(y)-1) = quadtree;

for i = 1:length(x)-1
    for j = 1:length(y)-1
        r(i,j) = quadtree(0, x(i), x(i+1), y(j), y(j+1), true, quadtree.empty, 'root');
    end
end


for i = 1:length(x)-1
    for j = 1:length(y)-1
        r(i,j).create_tree([1;1], 2, 'points')
        r(i,j).balance_tree(1);
        r(i,j).draw_tree()
        r(i,j).add_leafs()
    end
end

axis image

leafs = [r(:,:).unbalanced_leafs];

for leaf = leafs
    leaf.find_neighbors()
end

%%
hold on

nx = 128;
ny = 128;

tstep = 50000;

% nu = 0.0005;
tau = 0.6;
% tau = (6*dt*nu/dx^2 + 1)/2
omega = 1/tau;
% nu = (2*tau - 1)/6*dx^2/dt;

u_ini = 0.3;
v_ini = 0;

% density = 1.0;

drawstep = 1;

w1 = 4/9;
w2 = 1/9;
w3 = 1/36;

% f1 = gpuArray(w2*ones(nx,ny));

for i = 1:size(leafs,2)
    xx(i,:) = linspace(leafs(i).xmin, leafs(i).xmax, nx);
    yy(i,:) = linspace(leafs(i).ymin, leafs(i).ymax, ny);
    
    dx(i) = xx(i,2)-xx(i,1);
    dy(i) = dx(i);
    
    dt(i) = dx(i);
    c = dx./dt;
    
    f1 = w2*ones(nx+2,ny+2,i);
    f2 = w3*ones(nx+2,ny+2,i);
    f3 = w2*ones(nx+2,ny+2,i);
    f4 = w3*ones(nx+2,ny+2,i);
    f5 = w2*ones(nx+2,ny+2,i);
    f6 = w3*ones(nx+2,ny+2,i);
    f7 = w2*ones(nx+2,ny+2,i);
    f8 = w3*ones(nx+2,ny+2,i);
    f9 = w1*ones(nx+2,ny+2,i);
end

f_eq1 = f1;
f_eq2 = f2;
f_eq3 = f3;
f_eq4 = f4;
f_eq5 = f5;
f_eq6 = f6;
f_eq7 = f7;
f_eq8 = f8;
f_eq9 = f9;

%%

for t = 1:tstep
    
    for i = [1 2 4]
        
        % interior
        f1(2:nx, 2:ny, i) = f1(1:nx-1, 2:ny  , i);
        f2(2:nx, 2:ny, i) = f2(1:nx-1, 1:ny-1, i);
        f3(2:nx, 2:ny, i) = f3(2:nx  , 1:ny-1, i);
        f4(2:nx, 2:ny, i) = f4(3:nx+1, 1:ny-1, i);
        f5(2:nx, 2:ny, i) = f5(3:nx+1, 2:ny  , i);
        f6(2:nx, 2:ny, i) = f6(3:nx+1, 3:ny+1, i);
        f7(2:nx, 2:ny, i) = f7(2:nx  , 3:ny+1, i);
        f8(2:nx, 2:ny, i) = f8(1:nx-1, 3:ny+1, i);
        
        % Boundary
        % i=1, Bounceback
        if i == 1 || i == 3
            f1(2, 2:ny, i) = f5(2, 2:ny, i);
            f2(2, 2:ny, i) = f6(2, 2:ny, i);
            f8(2, 2:ny, i) = f4(2, 2:ny, i);
        end
        
        if i == 2 || i == 4
            % i=nx, Bounceback
            f4(nx, 2:ny, i) = f8(nx, 2:ny, i);
            f5(nx, 2:ny, i) = f1(nx, 2:ny, i);
            f6(nx, 2:ny, i) = f2(nx, 2:ny, i);
        end
        
        if i == 3 || i == 4
            % j=1, Bounceback
            f2(2:nx, 2, i) = f6(2:nx, 2, i);
            f3(2:nx, 2, i) = f7(2:nx, 2, i);
            f4(2:nx, 2, i) = f8(2:nx, 2, i);
        end
        
        if i == 1 || i == 2
            % j=ny, Know Velocity
            densityN = f9(2:nx, ny, i) + f1(2:nx, ny, i) + f5(2:nx, ny, i) + 2*(f3(2:nx, ny, i) + f4(2:nx, ny, i) + f2(2:nx, ny, i));
            
            f6(2:nx, ny, i) = f2(2:nx, ny, i) + 0.5*(f1(2:nx, ny, i) - f5(2:nx, ny, i)) - 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
            f7(2:nx, ny, i) = f3(2:nx, ny, i) + 2/3*densityN.*v_ini;
            f8(2:nx, ny, i) = f4(2:nx, ny, i) + 0.5*(f5(2:nx, ny, i) - f1(2:nx, ny, i)) + 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
        end
        
        density = f1(2:nx, 2:ny, i) + f2(2:nx, 2:ny, i) + f3(2:nx, 2:ny, i) + f4(2:nx, 2:ny, i) + f5(2:nx, 2:ny, i) + f6(2:nx, 2:ny, i) + f7(2:nx, 2:ny, i) + f8(2:nx, 2:ny, i) + f9(2:nx, 2:ny, i);
        u = (f1(2:nx, 2:ny, i) + f2(2:nx, 2:ny, i) + f8(2:nx, 2:ny, i) - (f4(2:nx, 2:ny, i) + f5(2:nx, 2:ny, i) + f6(2:nx, 2:ny, i)))./density;
        v = (f2(2:nx, 2:ny, i) + f3(2:nx, 2:ny, i) + f4(2:nx, 2:ny, i) - (f6(2:nx, 2:ny, i) + f7(2:nx, 2:ny, i) + f8(2:nx, 2:ny, i)))./density;
        
        f_eq1(2:nx, 2:ny, i) = w2*density.*(1 + 3*u/c(i)      + 9/2*u.^2/c(i)^2      - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq3(2:nx, 2:ny, i) = w2*density.*(1 + 3*v/c(i)      + 9/2*v.^2/c(i)^2      - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq5(2:nx, 2:ny, i) = w2*density.*(1 - 3*u/c(i)      + 9/2*u.^2/c(i)^2      - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq7(2:nx, 2:ny, i) = w2*density.*(1 - 3*v/c(i)      + 9/2*v.^2/c(i)^2      - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq2(2:nx, 2:ny, i) = w3*density.*(1 + 3*( u+v)/c(i) + 9/2*( u+v).^2/c(i)^2 - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq4(2:nx, 2:ny, i) = w3*density.*(1 + 3*(-u+v)/c(i) + 9/2*(-u+v).^2/c(i)^2 - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq6(2:nx, 2:ny, i) = w3*density.*(1 - 3*( u+v)/c(i) + 9/2*( u+v).^2/c(i)^2 - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq8(2:nx, 2:ny, i) = w3*density.*(1 + 3*( u-v)/c(i) + 9/2*( u-v).^2/c(i)^2 - 3/2*(u.^2 + v.^2)/c(i)^2);
        f_eq9(2:nx, 2:ny, i) = w1*density.*(1                                        - 3/2*(u.^2 + v.^2)/c(i)^2);
        
        f1(2:nx, 2:ny, i) = omega*f_eq1(2:nx, 2:ny, i) + (1-omega)*f1(2:nx, 2:ny, i);
        f2(2:nx, 2:ny, i) = omega*f_eq2(2:nx, 2:ny, i) + (1-omega)*f2(2:nx, 2:ny, i);
        f3(2:nx, 2:ny, i) = omega*f_eq3(2:nx, 2:ny, i) + (1-omega)*f3(2:nx, 2:ny, i);
        f4(2:nx, 2:ny, i) = omega*f_eq4(2:nx, 2:ny, i) + (1-omega)*f4(2:nx, 2:ny, i);
        f5(2:nx, 2:ny, i) = omega*f_eq5(2:nx, 2:ny, i) + (1-omega)*f5(2:nx, 2:ny, i);
        f6(2:nx, 2:ny, i) = omega*f_eq6(2:nx, 2:ny, i) + (1-omega)*f6(2:nx, 2:ny, i);
        f7(2:nx, 2:ny, i) = omega*f_eq7(2:nx, 2:ny, i) + (1-omega)*f7(2:nx, 2:ny, i);
        f8(2:nx, 2:ny, i) = omega*f_eq8(2:nx, 2:ny, i) + (1-omega)*f8(2:nx, 2:ny, i);
        f9(2:nx, 2:ny, i) = omega*f_eq9(2:nx, 2:ny, i) + (1-omega)*f9(2:nx, 2:ny, i);
        
        if mod(t, 250) == 0
            imagesc(xx(i, 2:nx), yy(i, 2:ny), sqrt(u'.^2 + v'.^2))
%             contour(xx(i, 2:nx), yy(i, 2:ny), sqrt(u'.^2 + v'.^2), 100);
            
            title(['step is ',sprintf("%f",t)])
            axis image
            
            colormap(jet(64))
            % colorbar
            
            drawnow
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%   i == 1   %%%%%%%%%%%%%%%%%%%%
    % f1(1   , 2:ny, 1) = f1(nx  , 2:ny, 1);
    
    % f2(1   , 2:ny, 1) = f2(nx  , 2:ny, 3);
    f2(2:nx, 1   , 1) = f2(2:nx, ny  , 3);
    % f2(1   , 1   , 1) = f2(nx  , ny  , 1);
    
    f3(2:nx, 1   , 1) = f3(2:nx, ny  , 3);
    
    f4(nx+1, 2:ny, 1) = f4(2   , 2:ny, 2);
    f4(2:nx, 1   , 1) = f4(2:nx, ny  , 3);
    f4(nx+1, 1   , 1) = f4(2   , ny  , 4);
    
    f5(nx+1, 2:ny, 1) = f5(2   , 2:ny, 2);
    
    f6(nx+1, 2:ny, 1) = f6(2   , 2:ny, 2);
    %     f6(2:nx, ny+1, 1) = f6(2:nx, 2   , 2);
    %     f6(nx+1, ny+1, 1) = f6(2   , 2   , 2);
    
    %     f7(2:nx, ny+1, 1) = f7(2:nx, 2   , 1);
    
    %     f8(1   , 2:ny, 1) = f8(nx  , 2:ny, 1);
    %     f8(2:nx, ny+1, 1) = f8(2:nx, 2   , 1);
    %     f8(1   , ny+1, 1) = f8(nx  , 2   , 1);
    
    
    %%%%%%%%%%%%%%%%%%%%   i == 2   %%%%%%%%%%%%%%%%%%%%
    f1(1   , 2:ny, 2) = f1(nx  , 2:ny, 1);
    
    f2(1   , 2:ny, 2) = f2(nx  , 2:ny, 1);
    f2(2:nx, 1   , 2) = f2(2:nx, ny  , 4);
    f2(1   , 1   , 2) = f2(nx  , ny  , 3);
    
    f3(2:nx, 1   , 2) = f3(2:nx, ny  , 4);
    
%     f4(nx+1, 2:ny, 2) = f4(2   , 2:ny, 4);
    f4(2:nx, 1   , 2) = f4(2:nx, ny  , 4);
%     f4(nx+1, 1   , 2) = f4(2   , ny  , 3);
    
    %     f5(nx    , 2:ny-1, 2) = f5(1     , 2:ny-1, 2);
    
    %     f6(nx    , 2:ny-1, 2) = f6(1     , 2:ny-1, 2);
    %     f6(2:nx-1, ny    , 2) = f6(2:nx-1, 2     , 1);
    %     f6(nx    , ny    , 2) = f6(2     , 2     , 3);
    
    %     f7(2:nx-1, ny    , 2) = f7(2:nx-1, 2     , 2);
    
    f8(1  , 2:ny, 2) = f8(nx   , 2:ny, 1);
    %     f8(2:nx-1, ny    , 2) = f8(2:nx-1, 2     , 1);
    %     f8(1  , ny+1, 2) = f8(nx   , 2   , 1);
    
    %%%%%%%%%%%%%%%%%%%%   i == 3   %%%%%%%%%%%%%%%%%%%%
    %     f1(1   , 2:ny, 3) = f1(nx  , 2:ny, 1);
    
    %     f2(1   , 2:ny, 3) = f2(nx  , 2:ny, 3);
    %     f2(2:nx, 1   , 3) = f2(2:nx, ny  , 1);
    %     f2(1   , 1   , 3) = f2(nx  , ny  , 1);
    
    %     f3(2:nx, 1   , 3) = f3(2:nx, ny  , 3);
    
    f4(nx+1, 2:ny, 3) = f4(2   , 2:ny, 4);
    %     f4(2:nx, 1   , 3) = f4(2:nx, ny  , 3);
    %     f4(nx+1, 1   , 3) = f4(2   , ny  , 4);
    
    f5(nx+1, 2:ny, 3) = f5(2   , 2:ny, 4);
    
    f6(nx+1, 2:ny, 3) = f6(2   , 2:ny, 4);
    f6(2:nx, ny+1, 3) = f6(2:nx, 2   , 1);
    f6(nx+1, ny+1, 3) = f6(2   , 2   , 2);
    
    f7(2:nx, ny+1, 3) = f7(2:nx, 2   , 1);
    
    %     f8(1   , 2:ny, 3) = f8(nx  , 2:ny, 1);
    f8(2:nx, ny+1, 3) = f8(2:nx, 2   , 1);
    %     f8(1   , ny+1, 3) = f8(nx  , 2   , 1);
    
    %%%%%%%%%%%%%%%%%%%%   i == 4   %%%%%%%%%%%%%%%%%%%%
    f1(1   , 2:ny, 4) = f1(nx  , 2:ny, 3);
    
    f2(1   , 2:ny, 4) = f2(nx  , 2:ny, 3);
    %     f2(2:nx, 1   , 4) = f2(2:nx, ny  , 1);
    %     f2(1   , 1   , 4) = f2(nx  , ny  , 1);
    
    %     f3(2:nx, 1   , 4) = f3(2:nx, ny  , 3);
    
    %     f4(nx+1, 2:ny, 4) = f4(2   , 2:ny, 2);
    %     f4(2:nx, 1   , 4) = f4(2:nx, ny  , 3);
    %     f4(nx+1, 1   , 4) = f4(2   , ny  , 4);
    
    %     f5(nx+1, 2:ny, 4) = f5(2   , 2:ny, 2);
    
    %     f6(nx+1, 2:ny, 4) = f6(2   , 2:ny, 2);
    f6(2:nx, ny+1, 4) = f6(2:nx, 2   , 2);
    %     f6(nx+1, ny+1, 4) = f6(2   , 2   , 2);
    
    f7(2:nx, ny+1, 4) = f7(2:nx, 2   , 2);
    
    f8(1   , 2:ny, 4) = f8(nx  , 2:ny, 3);
    f8(2:nx, ny+1, 4) = f8(2:nx, 2   , 2);
    f8(1   , ny+1, 4) = f8(nx  , 2   , 1);
    
end

hold off

































