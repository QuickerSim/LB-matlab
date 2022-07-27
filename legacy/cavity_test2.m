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
        r(i,j).create_tree([1; 1], 3, 'points')
        r(i,j).balance_tree(1);
%         r(i,j).draw_tree()
        r(i,j).add_leafs()
    end
end

leafs = [r(:,:).unbalanced_leafs];

for leaf = leafs
    leaf.find_neighbors()
end

disp('done')

%%

% leafs(1).initLB_TestStreaming(128, 128, 0.6, 0.3, 0);
% for i = 1:size(leafs,2)
%     leafs(i).initLB(128, 128, 0.6, 0.3, 0);
% end
%
% % for i = 1:size(leafs,2)
% %     leafs(i).initLB(256, 256, 0.7, 0.1, 0);
% % end
%
%
% for t = 1:5000
%     disp(t)
%
%     for i = 1:size(leafs,2)
%         leafs(i).Collision();
%     end
%
%     for i = 1:size(leafs,2)
%         leafs(i).GhostInfo();
%     end
%
%     for i = 1:size(leafs,2)
%         leafs(i).Streaming();
%     end
%
%     if mod(t,5000) == 0
%         figure(44)
%         clf
%         hold on
%         for i = 1:size(leafs,2)
%             imagesc(leafs(i).x, leafs(i).y, sqrt(leafs(i).u.^2 + leafs(i).v.^2)');
%             %             contourf(leafs(i).x, leafs(i).y, sqrt(leafs(i).u.^2 + leafs(i).v.^2)', 0:0.0005:0.001);
%             leafs(i).draw_tree()
%         end
%         axis image
%         caxis([0 0.3])
%         colormap(jet(5555))
%         drawnow
%         hold off
%     end
%
% end

%%

% figure(44)
%
% clf
%
% hold on
%
% for i = 1:size(leafs,2)
%     contourf(leafs(i).x, leafs(i).y, sqrt(leafs(i).u.^2 + leafs(i).v.^2)',0:0.01:0.25)
%     leafs(i).draw_tree()
% end
%
% axis image
%
% hold off

%%

sor = sortrows([1:size(leafs,2) ;leafs.deepness]',2);

leafs = leafs(sor(:,1));

leafs.deepness

clc

leafs(1).initLB_TestStreaming(64, 64, 0.7, 0.3, 0);
for i = 2:size(leafs,2)
    leafs(i).initLB(64, 64, 0.7, 0.3, 0);
end

% for i = 1:size(leafs,2)
%     leafs(i).initLB(64, 64, 0.7, 0.3, 0);
% end

for t = 1:10000
    disp(t)
    
    lbdo(leafs,1)
    
    if mod(t,5) == 0
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
        caxis([0 0.001])
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

for i = ind
    leafs(i).Streaming();
end

lbdo(leafs, deepness+1)

for i = ind
    leafs(i).InterpolateUp();
end

end






















