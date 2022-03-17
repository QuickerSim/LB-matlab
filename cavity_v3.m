clear all
clc

x = 0:1:1;
y = 0:1:1;

r(length(x)-1, length(y)-1) = quadtreelb;

for i = 1:length(x)-1
    for j = 1:length(y)-1
        r(i,j) = quadtreelb(0, x(i), x(i+1), y(j), y(j+1), true, quadtreelb.empty, 'root');
    end
end

for i = 1:length(x)-1
    for j = 1:length(y)-1
        r(i,j).create_tree([0 1; 1 1], 1, 'points')
        r(i,j).balance_tree(1);
        r(i,j).add_leafs()
    end
end

leafs = [r(:,:).unbalanced_leafs];

for leaf = leafs
    leaf.find_neighbors()
end

disp('done')

%%

sor = sortrows([1:size(leafs,2) ;leafs.deepness]',2);

leafs = leafs(sor(:,1));

leafs.deepness

clc

n = 256;

for i = 1:size(leafs,2)
    leafs(i).initLB(n, n, 0.7, 0.1, 0.01);
end

% leafs(36).initLB_TestStreaming(n, n, 0.7, 0.3, 0);

steps = 50000;

for t = 1:steps
    disp(t)
    
    lbdo(leafs, min([leafs.deepness]))
    
    if mod(t,50) == 0
        figure(44)
        clf
        hold on
        for i = 1:size(leafs,2)
            imagesc(leafs(i).x, leafs(i).y, sqrt(leafs(i).u.^2 + leafs(i).v.^2)');
            % contourf(leafs(i).x, leafs(i).y, sqrt(leafs(i).u.^2 + leafs(i).v.^2)', 0:0.0005:0.001);
            leafs(i).draw_tree()
        end
        axis image
        colormap(jet(5555))
%         caxis([0 0.05])
        drawnow
        hold off
    end
    
end


function lbdo(leafs, deepness)

ind = find([leafs.deepness] == deepness);

if isempty(ind)
    return
end

% main loop
for i = ind
    leafs(i).InterpolateDown();
end

for i = ind
    leafs(i).Collision();
end

lbdo(leafs, deepness+1)

for i = ind
    leafs(i).GhostInfo();
end

lbdo(leafs, deepness+1)

for i = ind
    leafs(i).Streaming();
end

for i = ind
    leafs(i).InterpolateUp();
end

end






















