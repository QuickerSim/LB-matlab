tic; hold on; clc; clear; nx=100; ny=100; tstep=40000; 
alpha=0.01; omega=1.0; u_ini=0.1; v_ini=0; Re=u_ini*nx/alpha; 
w1=4/9; w2=1/9; w3=1/36; 
f=ones(nx,ny,9); 
f_eq=f; 
density=2.7;

for ii=1:2

    % Propegate (This part of code [propegate] is always constant for all LBM
    % problems.)
    f(:,:,4)=f([2:nx 1],[ny 1:ny-1],4); f(:,:,3)=f(:,[ny 1:ny-1],3);
    f(:,:,2)=f([nx 1:nx-1],[ny 1:ny-1],2); f(:,:,5)=f([2:nx 1],:,5);
    f(:,:,1)=f([nx 1:nx-1],:,1); f(:,:,6)=f([2:nx 1],[2:ny 1],6);
    f(:,:,7)=f(:,[2:ny 1],7); f(:,:,8)=f([nx 1:nx-1],[2:ny 1],8);
    % Boundary Conditions
    %At i=1, Bounceback
    f(1,:,1)=f(1,:,5); f(1,:,2)=f(1,:,6); f(1,:,8)=f(1,:,4);
    %At i=nx, Bounceback
    f(nx,:,4)=f(nx,:,8); f(nx,:,5)=f(nx,:,1); f(nx,:,6)=f(nx,:,2);
    %At j=1, Bounceback
    f(:,1,2)=f(:,1,6); f(:,1,3)=f(:,1,7); f(:,1,4)=f(:,1,8);
    %At j=ny, Know Velocity
    densityN=f(:,ny,9)+f(:,ny,1)+f(:,ny,5)+2*(f(:,ny,3)+f(:,ny,4)+f(:,ny,2));
    f(:,ny,7)=f(:,ny,3);
    f(:,ny,6)=f(:,ny,2)+0.5*(f(:,ny,1)-f(:,ny,5))-u_ini.*densityN/2;
    f(:,ny,8)=f(:,ny,4)+0.5*(f(:,ny,5)-f(:,ny,1))+u_ini.*densityN/2;
    density=sum(f,3);
    u=(sum(f(:,:,[1 2 8]),3)-sum(f(:,:,[4 5 6]),3))./density;
    v=(sum(f(:,:,[2 3 4]),3)-sum(f(:,:,[6 7 8]),3))./density;
    f_eq(:,:,1)=w2*density.*(1+3*u+9/2*u.^2-3/2*(u.^2+v.^2));
    f_eq(:,:,3)=w2*density.*(1+3*v+9/2*v.^2-3/2*(u.^2+v.^2));
    f_eq(:,:,5)=w2*density.*(1-3*u+9/2*u.^2-3/2*(u.^2+v.^2));
    f_eq(:,:,7)=w2*density.*(1-3*v+9/2*v.^2-3/2*(u.^2+v.^2));
    f_eq(:,:,2)=w3*density.*(1+3*(u+v)+9/2*(u+v).^2-3/2*(u.^2+v.^2));
    f_eq(:,:,4)=w3*density.*(1+3*(-u+v)+9/2*(-u+v).^2-3/2*(u.^2+v.^2));
    f_eq(:,:,6)=w3*density.*(1-3*(u+v)+9/2*(u+v).^2-3/2*(u.^2+v.^2));
    f_eq(:,:,8)=w3*density.*(1+3*(u-v)+9/2*(u-v).^2-3/2*(u.^2+v.^2));
    f_eq(:,:,9)=w1*density.*(1-3/2*(u.^2+v.^2));
    f=omega*f_eq+(1-omega)*f;    
    imagesc(u')
    drawnow
end

contour(u',100); toc;