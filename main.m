clear all
clc
clf

% mesh
[p,e,t] = importMeshGmsh('cyl.msh');
boundary = unique(e([1 2],:));
pbound = p(:, boundary);

% quadtree
root = quadtree(0, 0, 5, 0, 1, true, quadtree.empty, 'root');
root.create_tree(pbound, 5)
root.draw_tree()
axis image

%%

root.balance_tree(2)
disp('tree balanced')

clf
root.draw_tree()
axis image
disp('tree ploted')

%% forest

clear r
clf
% x = linspace(0,5,11)
% y = linspace(0,1,3)

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
%         r(i,j).balance_tree(1);
        r(i,j).draw_tree()
        r(i,j).add_leafs()          
    end
end

axis image

disp('done')



leafs = [r(:,:).unbalanced_leafs];
% test = find([leafs.deepness] == 2);
% for i = leafs(test)
%    i.draw_tree('red') 
% end

%%


for i = 1:length(x)-1
    for j = 1:length(y)-1
        
    end
end





























