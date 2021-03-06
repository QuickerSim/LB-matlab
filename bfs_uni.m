clc;
clear all;
% close all;

a = 2;
b = 1;

nx = 200;
ny = 200;

x = linspace(0,a,nx);
y = linspace(0,b,ny);

dx = x(2)-x(1);
dy = dx;

nx = length(x)
ny = length(y)

tstep = 100000;
dt = dx;

c = dx/dt

nu = 0.0001;
tau = (6*dt*nu/dx^2 + 1)/2
omega = 1/tau;
nu = (2*tau - 1)/6*dx^2/dt

u_ini = 0;
v_ini = 0.1;

alpha = 0.1;

Re = u_ini*nx/alpha;

density = 1.0;

drawstep = 1000;


%%
w1 = 4/9;
w2 = 1/9;
w3 = 1/36;

f1 = w2*ones(nx,ny);
f2 = w3*ones(nx,ny);
f3 = w2*ones(nx,ny);
f4 = w3*ones(nx,ny);
f5 = w2*ones(nx,ny);
f6 = w3*ones(nx,ny);
f7 = w2*ones(nx,ny);
f8 = w3*ones(nx,ny);
f9 = w1*ones(nx,ny);

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

for ii = 1:tstep
    tic
    % interior
    f1(2:nx-1, 2:ny-1) = f1([1:nx-2], [2:ny-1]);
    f2(2:nx-1, 2:ny-1) = f2([1:nx-2], [1:ny-2]);
    f3(2:nx-1, 2:ny-1) = f3([2:nx-1], [1:ny-2]);
    f4(2:nx-1, 2:ny-1) = f4([3:nx]  , [1:ny-2]);
    f5(2:nx-1, 2:ny-1) = f5([3:nx]  , [2:ny-1]);
    f6(2:nx-1, 2:ny-1) = f6([3:nx]  , [3:ny]  );
    f7(2:nx-1, 2:ny-1) = f7([2:nx-1], [3:ny]  );
    f8(2:nx-1, 2:ny-1) = f8([1:nx-2], [3:ny]  );
    
    f1(1:nx/8, 1:ny/4) = w2;
    f2(1:nx/8, 1:ny/4) = w3;
    f3(1:nx/8, 1:ny/4) = w2;
    f4(1:nx/8, 1:ny/4) = w3;
    f5(1:nx/8, 1:ny/4) = w2;
    f6(1:nx/8, 1:ny/4) = w3;
    f7(1:nx/8, 1:ny/4) = w2;
    f8(1:nx/8, 1:ny/4) = w3;
    f9(1:nx/8, 1:ny/4) = w1;
    
    % Boundary
    % i=1, Bounceback - left
    densityN = f9(1,ny/4+1:end-1) + f3(1,ny/4+1:end-1) + f7(1,ny/4+1:end-1) + 2*(f4(1,ny/4+1:end-1) + f5(1,ny/4+1:end-1) + f6(1,ny/4+1:end-1));
    
    f2(1,ny/4+1:end-1) = f6(1,ny/4+1:end-1) + 0.5*(f7(1,ny/4+1:end-1) - f3(1,ny/4+1:end-1)) - 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    f1(1,ny/4+1:end-1) = f5(1,ny/4+1:end-1) + 2/3*densityN.*v_ini;
    f8(1,ny/4+1:end-1) = f4(1,ny/4+1:end-1) + 0.5*(f3(1,ny/4+1:end-1) - f7(1,ny/4+1:end-1)) + 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    
    % i=nx, Bounceback - right
    % f4(nx,:) = f8(nx,:);
    % f5(nx,:) = f1(nx,:);
    % f6(nx,:) = f2(nx,:);
    densityN = f9(nx,:) + f3(nx,:) + f7(nx,:) + 2*(f4(nx,:) + f5(nx,:) + f6(nx,:));
    densityN = 0;
    
    f2(nx,:) = f6(nx,:) + 0.5*(f7(nx,:) - f3(nx,:)) - 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    f1(nx,:) = f5(nx,:) + 2/3*densityN.*v_ini;
    f8(nx,:) = f4(nx,:) + 0.5*(f3(nx,:) - f7(nx,:)) + 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    
    % bfs, Bounceback - step - right
    f4(nx/8,1:ny/4) = f8(nx/8,1:ny/4);
    f5(nx/8,1:ny/4) = f1(nx/8,1:ny/4);
    f6(nx/8,1:ny/4) = f2(nx/8,1:ny/4);
    
    % bfs, Bounceback - step - top
    f2(1:nx/8,ny/4) = f6(1:nx/8,ny/4);
    f3(1:nx/8,ny/4) = f7(1:nx/8,ny/4);
    f4(1:nx/8,ny/4) = f8(1:nx/8,ny/4);
    
    % j=1, Bounceback - bottom
    f2(:,1) = f6(:,1);
    f3(:,1) = f7(:,1);
    f4(:,1) = f8(:,1);
    
    % j=ny, Bounceback - top
    f6(:,ny) = f2(:,ny);
    f7(:,ny) = f3(:,ny);
    f8(:,ny) = f4(:,ny);
    
    density = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9;
    u = (f1 + f2 + f8 - (f4 + f5 + f6))./density;
    v = (f2 + f3 + f4 - (f6 + f7 + f8))./density;
    
    f_eq1 = w2*density.*(1 + 3*u/c      + 9/2*u.^2/c^2      - 3/2*(u.^2 + v.^2)/c^2);
    f_eq3 = w2*density.*(1 + 3*v/c      + 9/2*v.^2/c^2      - 3/2*(u.^2 + v.^2)/c^2);
    f_eq5 = w2*density.*(1 - 3*u/c      + 9/2*u.^2/c^2      - 3/2*(u.^2 + v.^2)/c^2);
    f_eq7 = w2*density.*(1 - 3*v/c      + 9/2*v.^2/c^2      - 3/2*(u.^2 + v.^2)/c^2);
    f_eq2 = w3*density.*(1 + 3*( u+v)/c + 9/2*( u+v).^2/c^2 - 3/2*(u.^2 + v.^2)/c^2);
    f_eq4 = w3*density.*(1 + 3*(-u+v)/c + 9/2*(-u+v).^2/c^2 - 3/2*(u.^2 + v.^2)/c^2);
    f_eq6 = w3*density.*(1 - 3*( u+v)/c + 9/2*( u+v).^2/c^2 - 3/2*(u.^2 + v.^2)/c^2);
    f_eq8 = w3*density.*(1 + 3*( u-v)/c + 9/2*( u-v).^2/c^2 - 3/2*(u.^2 + v.^2)/c^2);
    f_eq9 = w1*density.*(1                                  - 3/2*(u.^2 + v.^2)/c^2);
    
    f1 = omega*f_eq1 + (1-omega)*f1;
    f2 = omega*f_eq2 + (1-omega)*f2;
    f3 = omega*f_eq3 + (1-omega)*f3;
    f4 = omega*f_eq4 + (1-omega)*f4;
    f5 = omega*f_eq5 + (1-omega)*f5;
    f6 = omega*f_eq6 + (1-omega)*f6;
    f7 = omega*f_eq7 + (1-omega)*f7;
    f8 = omega*f_eq8 + (1-omega)*f8;
    f9 = omega*f_eq9 + (1-omega)*f9;
    toc
    if mod(ii,drawstep) == 0
        %         surf(1:nx,ny:-1:1,sqrt(u'.^2 + v'.^2),'EdgeColor','none')
        %         view([210 30])
        %                 imagesc(x,ny:-1:1,sqrt(u'.^2 + v'.^2))
%         contourf(x,y,sqrt(u'.^2 + v'.^2),20);
        contour(x,y,sqrt(u'.^2 + v'.^2),100);
        %         contourf(x,y,density');
        
        title(['step is ',sprintf("%f",ii)])
        axis image
        colormap(jet(64))
        colorbar
        caxis([0 0.1])
        drawnow
        
        
    end
end













