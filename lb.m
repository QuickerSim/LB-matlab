clc;
clear all;
% close all;


a = 2;
b = 1;

nx = 200;
ny = 100;

x = linspace(0,a,nx);
y = linspace(0,b,ny);

dx = x(2)-x(1)
dy = dx;

nx = length(x)
ny = length(y)

tstep = 10000;
dt = dx;

c = dx/dt

nu = 0.001;
% tau = 1;
tau = (6*dt*nu/dx^2 + 1)/2
omega = 1/tau
nu = (2*tau - 1)/6*dx^2/dt

u_ini = 0.3;
v_ini = 0;

alpha = 0.1;

Re = u_ini*nx/alpha;

density = 1.0;

drawstep = 100;


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
    
    % interior
    f1(1:nx, 1:ny) = f1([nx 1:nx-1], [1:ny]     );
    f2(1:nx, 1:ny) = f2([nx 1:nx-1], [ny 1:ny-1]);
    f3(1:nx, 1:ny) = f3([1:nx]     , [ny 1:ny-1]);
    f4(1:nx, 1:ny) = f4([2:nx 1]   , [ny 1:ny-1]);
    f5(1:nx, 1:ny) = f5([2:nx 1]   , [1:ny]     );
    f6(1:nx, 1:ny) = f6([2:nx 1]   , [2:ny 1]   );
    f7(1:nx, 1:ny) = f7([1:nx]     , [2:ny 1]   );
    f8(1:nx, 1:ny) = f8([nx 1:nx-1], [2:ny 1]   );    
    
    % Boundary
    
    % i=1, Bounceback
    f1(1,:) = f5(1,:);
    f2(1,:) = f6(1,:);
    f8(1,:) = f4(1,:);
    
    % i=nx, Bounceback
    f4(nx,:) = f8(nx,:);
    f5(nx,:) = f1(nx,:);
    f6(nx,:) = f2(nx,:);
    
    % bfs? step - right
    f4(nx/4,1:ny/2) = f8(nx/4,1:ny/2);
    f5(nx/4,1:ny/2) = f1(nx/4,1:ny/2);
    f6(nx/4,1:ny/2) = f2(nx/4,1:ny/2);
    
    % bfs? step - top
    f2(1:nx/4,ny/2) = f6(1:nx/4,ny/2);
    f3(1:nx/4,ny/2) = f7(1:nx/4,ny/2);
    f4(1:nx/4,ny/2) = f8(1:nx/4,ny/2);
    
    % j=1, Bounceback
    f2(:,1) = f6(:,1);
    f3(:,1) = f7(:,1);
    f4(:,1) = f8(:,1);
    
    % j=ny, Know Velocity
    densityN = f9(:,ny) + f1(:,ny) + f5(:,ny) + 2*(f3(:,ny) + f4(:,ny) + f2(:,ny));
   
    f6(:,ny) = f2(:,ny) + 0.5*(f1(:,ny) - f5(:,ny)) - 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    f7(:,ny) = f3(:,ny) + 2/3*densityN.*v_ini;
    f8(:,ny) = f4(:,ny) + 0.5*(f5(:,ny) - f1(:,ny)) + 0.5*u_ini.*densityN + 1/6*densityN*v_ini;
    
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
    
    if mod(ii,drawstep) == 0
%         surf(1:nx,ny:-1:1,sqrt(u'.^2 + v'.^2),'EdgeColor','none')
%         view([210 30])
%         imagesc(x,y,sqrt(u'.^2 + v'.^2))
        contourf(x,y,sqrt(u'.^2 + v'.^2),20);
        
        title(['step is ',sprintf("%f",ii)])
        axis image
        drawnow
        
        colormap(jet(64))
        colorbar
    end
end
%%
contourf(sqrt(u'.^2 + v'.^2),20);












