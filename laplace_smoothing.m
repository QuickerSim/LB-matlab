clear all
clc

%%

x = 1:100;
y = 1:100;

[X, Y] = meshgrid(x,y);

init = @(x,y) x.^2 + y.^2 + 200000*cos(x).*sin(y);

A = init(X,Y);

surf(A)

%%
% for i = 1:10000
%
%     A(1,1)     = (A(1,1) + A(1,2) + A(2,1) + A(2,2))/4;
%     A(1,100)   = (A(1,100) + A(1,99) + A(2,100) + A(2,99))/4;
%     A(100,1)   = (A(100,1) + A(100,2) + A(99,1) + A(99,2))/4;
%     A(100,100) = (A(100,100) + A(100,99) + A(99,100) + A(99,99))/4;
%
%     A(2:99,1) = (A(1:98,1) + A(2:99,1) + A(3:100,1) + ...
%         A(1:98,2) + A(2:99,2) + A(3:100,2))/6;
%
%     A(2:99,100) = (A(1:98,100) + A(2:99,100) + A(3:100,100) + ...
%         A(1:98,99) + A(2:99,99) + A(3:100,99))/6;
%
%     A(1,2:99) = (A(1,1:98) + A(1,2:99) + A(1,3:100) + ...
%         A(2,1:98) + A(2,2:99) + A(2,3:100))/6;
%
%     A(100,2:99) = (A(100,1:98) + A(100,2:99) + A(100,3:100) + ...
%         A(99,1:98) + A(99,2:99) + A(99,3:100))/6;
%
%     A(2:99,2:99) = (A(1:98,1:98) + A(2:99,1:98) + A(3:100,1:98) + ...
%         A(1:98,2:99) + A(2:99,2:99) + A(3:100,2:99) + ...
%         A(1:98,3:100) + A(2:99,3:100) + A(3:100,3:100))/9;
%
%     if mod(i, 1) == 0
%         disp(i)
%         surf(A)
%         drawnow
%         %         pause(0.2)
%     end
% end
%% LBM from some matlab central

% interior
% f1(1:nx, 1:ny) = f1([nx 1:nx-1], [1:ny]     );
% f2(1:nx, 1:ny) = f2([nx 1:nx-1], [ny 1:ny-1]);
% f3(1:nx, 1:ny) = f3([1:nx]     , [ny 1:ny-1]);
% f4(1:nx, 1:ny) = f4([2:nx 1]   , [ny 1:ny-1]); 
% f5(1:nx, 1:ny) = f5([2:nx 1]   , [1:ny]     ); 
% f6(1:nx, 1:ny) = f6([2:nx 1]   , [2:ny 1]   );
% f7(1:nx, 1:ny) = f7([1:nx]     , [2:ny 1]   ); 
% f8(1:nx, 1:ny) = f8([nx 1:nx-1], [2:ny 1]   );


%%
Amat = sparse(size(x,2)*size(x,2),size(x,2)*size(x,2));

%top?
for i=1:size(x,2)
    if mod(i,size(x,2))~=1
        Amat(i-1,i) = 1;
        Amat(i-1+size(x,2),i) = 1;
    end
    
    Amat(i,i) = 1;
    Amat(i+size(x,2),i) = 1;
    
    if mod(i,size(x,2))~=0
        Amat(i+1,i) = 1;
        Amat(i+1+size(x,2),i) = 1;
    end
end

%bottom?
for i=size(x,2)*size(x,2)-size(x,2):size(x,2)*size(x,2)
    if mod(i,size(x,2))~=1
        Amat(i-1,i) = 1;
        Amat(i-1-size(x,2),i) = 1;
    end
    
    Amat(i,i) = 1;
    Amat(i-size(x,2),i) = 1;
    
    if mod(i,size(x,2))~=0
        Amat(i+1,i) = 1;
        Amat(i+1-size(x,2),i) = 1;
    end
end

%interior with left and right?
for i=size(x,2)+1:size(x,2)*size(x,2)-size(x,2)-1
    if mod(i,size(x,2))~=1
        Amat(i-1-size(x,2),i) = 1;
        Amat(i-1,i) = 1;
        Amat(i-1+size(x,2),i) = 1;
    end
    
    if mod(i,size(x,2))~=0
        Amat(i+1-size(x,2),i) = 1;
        Amat(i+1,i) = 1;
        Amat(i+1+size(x,2),i) = 1;
    end
    
    Amat(i-size(x,2),i) = 1;
    Amat(i,i) = 1;
    Amat(i+size(x,2),i) = 1;
end



%%
% Amat = Amat';
A = reshape(A,[],1);
%%
for i = size(x,2)*size(x,2):-1:1
    mian(i) = sum(Amat(i,:),2);
end
%%
for i = 1:1000
    A = Amat*A;
    A = A./mian';
    
    if mod(i,1) == 0
        disp(i)
        surf(reshape(A,size(x,2),size(x,2)));
        drawnow
        %         pause(0.2)
    end
end

















