clear all
close all
clc

fID = fopen('mesh/elemlist.lsb');
leafs = fread(fID, [2,Inf], 'int64');
leafs = leafs';
fclose(fID);

% fID = fopen('mesh/bnd.lsb');
% bnd_leafs = fread(fID, [26,size(leafs,1)], 'int64');
% bnd_leafs = bnd_leafs;
% fclose(fID);


elemnums = cumsum(8.^[0:7]);

L = 0;

for j = 1:size(leafs,1)
    if leafs(j,1) == 0
        L = 0;
    elseif (leafs(j,1) >= 1 && leafs(j,1) <= 8)
        L = 1;
    elseif (leafs(j,1) >= 9 && leafs(j,1) <= 72)
        L = 2;
    elseif (leafs(j,1) >= 73 && leafs(j,1) <= 584)
        L = 3;
    elseif (leafs(j,1) >= 585 && leafs(j,1) <= 4680)
        L = 4;
    elseif (leafs(j,1) >= 4681 && leafs(j,1) <= 37448)
        L = 5;
    elseif (leafs(j,1) >= 37449 && leafs(j,1) <= 299593)
        L = 6;
    end
    
    x(j,4) = L;
    
    for k = 1:3
        x(j,k) = bin2dec(join(string(fliplr(bitget(leafs(j,1)-elemnums(L), k:3:25))),''));
    end
end

% bnd_leafs
% ttt = x(find(leafs(:,2)==8),:)
% sum(bnd_leafs)

%%

bcleafs = x(find(leafs(:,2) == 8),:)

%%

x(7,:)
x(7,2) = x(7,2)-1;
x(7,:)

%%

tmp = dec2bin(x(7,3:-1:1))

tmp2 = reshape(tmp',1,[])

southleaf = bin2dec(tmp2)+elemnums(x(7,4))

%%

find(leafs(:,1) == southleaf)

find(leafs(:,1) == floor((109 - 1)/8))

leafs(5,:)
x(5,:)






























