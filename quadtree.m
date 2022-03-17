classdef quadtree < handle
    properties
        parent
        
        nw
        ne
        se
        sw
        
        nneighbor
        eneighbor
        sneighbor
        wneighbor
        
        side
        
        marked_for_ref
        
        deepness
        
%         unbalanced_leafs(1,:) quadtree
        unbalanced_leafs = [];
        
        xmin
        xmax
        ymin
        ymax
    end
    
    methods
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function qt = quadtree(deep, left, right, bottom, top, mark, parent, side)
            if nargin == 0
                qt.deepness = -1;
                
                qt.xmin = -1;
                qt.xmax = -1;
                qt.ymin = -1;
                qt.ymax = -1;
                
                qt.nw = quadtree.empty;
                qt.ne = quadtree.empty;
                qt.sw = quadtree.empty;
                qt.se = quadtree.empty;
                
                qt.marked_for_ref = -1;
                qt.parent = quadtree.empty;
                qt.side = "";
            else
                
                qt.deepness = deep;
                
                qt.xmin = left;
                qt.xmax = right;
                qt.ymin = bottom;
                qt.ymax = top;
                
                qt.nw = quadtree.empty;
                qt.ne = quadtree.empty;
                qt.sw = quadtree.empty;
                qt.se = quadtree.empty;
                
                qt.marked_for_ref = mark;
                qt.parent = parent;
                qt.side = side;
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function create_tree(qt, points, deep, crossing_method)
            if qt.deepness >= deep
                return;
            end
            
            if nargin <= 3
                crossing_method = 'points';
            end
            
            if strcmpi(crossing_method, 'lines')
                qt.is_crossed_poly(points);
            elseif strcmpi(crossing_method, 'points')
                qt.is_crossed(points);
            else
                error('### crossing method ###');
            end
            
            if qt.marked_for_ref == 1
                qt.create_childrens();
                qt.nw.create_tree(points, deep, crossing_method);
                qt.ne.create_tree(points, deep, crossing_method);
                qt.se.create_tree(points, deep, crossing_method);
                qt.sw.create_tree(points, deep, crossing_method);
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function output = is_crossed(qt, points)
            x = find(points(1,:) >= qt.xmin & points(1,:) <= qt.xmax);
            
            if ~isempty(x)
                if (any(points(2,x) >= qt.ymin & points(2,x) <= qt.ymax))
                    qt.marked_for_ref = 1;
                    output = 1;
                    return;
                end
            end
            
            output = 0;
            return;
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function output = is_crossed_poly(qt, points)
            polysquare = polyshape([qt.xmin qt.xmax qt.xmax qt.xmin], [qt.ymax qt.ymax qt.ymin qt.ymin]);
            
            in = intersect(polysquare, points');
            
            if ~isempty(in)
                qt.marked_for_ref = 1;
                output = 1;
                return;
            end
            
            output = 0;
            return;
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function create_childrens(qt)
            yhalf = (qt.ymax - qt.ymin)/2;
            xhalf = (qt.xmax - qt.xmin)/2;
            
            qt.nw = quadtree(qt.deepness + 1, qt.xmin, qt.xmin + xhalf, qt.ymin + yhalf, qt.ymax, 0, qt, 'nw');
            qt.ne = quadtree(qt.deepness + 1, qt.xmin + xhalf, qt.xmax, qt.ymin + yhalf, qt.ymax, 0, qt, 'ne');
            qt.sw = quadtree(qt.deepness + 1, qt.xmin, qt.xmin + xhalf, qt.ymin, qt.ymin + yhalf, 0, qt, 'sw');
            qt.se = quadtree(qt.deepness + 1, qt.xmin + xhalf, qt.xmax, qt.ymin, qt.ymin + yhalf, 0, qt, 'se');
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function draw_tree(qt, color)
            if nargin == 1
                color = 'black';
            end
            
            if ~isempty(qt.nw)
                draw_tree(qt.nw, color)
                draw_tree(qt.ne, color)
                draw_tree(qt.sw, color)
                draw_tree(qt.se, color)
            end
            
            line([qt.xmin qt.xmax qt.xmax qt.xmin qt.xmin], [qt.ymax qt.ymax qt.ymin qt.ymin qt.ymax], 'Color', color)
            
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function add_leafs(qt, root)
            if strcmp(qt.side, 'root')
                root = qt;
            end
            
            if ~isempty(qt.nw)
                qt.nw.add_leafs(root)
                qt.ne.add_leafs(root)
                qt.sw.add_leafs(root)
                qt.se.add_leafs(root)
            else
                root.unbalanced_leafs(end+1) = qt;
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function balance_tree(qt, deep_diff)
            if ~strcmp(qt.side, 'root')
                error('only root node can balance tree')
            end
            
            qt.add_leafs();
            
            while ~isempty(qt.unbalanced_leafs)
                
                leaf = qt.unbalanced_leafs(1);
                
                if strcmp(leaf.side, 'ne')
                    leaf_nn = leaf.north_neighbor();
                    if ~isempty(leaf_nn)
                        if leaf.deepness - leaf_nn.deepness > deep_diff
                            leaf_nn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_nn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.se;
                        end
                    end
                    
                    leaf_en = leaf.east_neighbor();
                    if ~isempty(leaf_en)
                        if leaf.deepness - leaf_en.deepness > deep_diff
                            leaf_en.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_en.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_en.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_en.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_en.se;
                        end
                    end
                elseif strcmp(leaf.side, 'nw')
                    leaf_nn = leaf.north_neighbor();
                    if ~isempty(leaf_nn)
                        if leaf.deepness - leaf_nn.deepness > deep_diff
                            leaf_nn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_nn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_nn.se;
                        end
                    end
                    
                    leaf_wn = leaf.west_neighbor();
                    if ~isempty(leaf_wn)
                        if leaf.deepness - leaf_wn.deepness > deep_diff
                            leaf_wn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_wn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.se;
                        end
                    end
                elseif strcmp(leaf.side, 'se')
                    leaf_sn = leaf.south_neighbor();
                    if ~isempty(leaf_sn)
                        if leaf.deepness - leaf_sn.deepness > deep_diff
                            leaf_sn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_sn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.se;
                        end
                    end
                    
                    leaf_en = leaf.east_neighbor();
                    if ~isempty(leaf_en)
                        if leaf.deepness - leaf_en.deepness > deep_diff
                            leaf_en.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_en.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_en.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_en.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_en.se;
                        end
                    end
                elseif strcmp(leaf.side, 'sw')
                    leaf_sn = leaf.south_neighbor();
                    if ~isempty(leaf_sn)
                        if leaf.deepness - leaf_sn.deepness > deep_diff
                            leaf_sn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_sn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_sn.se;
                        end
                    end
                    
                    leaf_wn = leaf.west_neighbor();
                    if ~isempty(leaf_wn)
                        if leaf.deepness - leaf_wn.deepness > deep_diff
                            leaf_wn.create_childrens();
                            qt.unbalanced_leafs(end + 1) = leaf_wn.nw;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.ne;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.sw;
                            qt.unbalanced_leafs(end + 1) = leaf_wn.se;
                        end
                    end
                end
                
                qt.unbalanced_leafs(1) = [];
            end
            
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function node = north_neighbor(qt)
            if strcmp(qt.side, 'root')
                node = quadtree.empty;
                return
            else
                par = qt.parent;
                if strcmp(qt.side, 'sw')
                    node = qt.parent.nw;
                    return
                elseif strcmp(qt.side, 'se')
                    node = qt.parent.ne;
                    return
                end
                mu = north_neighbor(qt.parent);
                if isempty(mu) || isempty(mu.nw)
                    node = mu;
                else
                    if strcmp(qt.side, 'nw')
                        node = mu.sw;
                        return
                    elseif strcmp(qt.side, 'ne')
                        node = mu.se;
                        return
                    end
                end
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function node = south_neighbor(qt)
            if strcmp(qt.side, 'root')
                node = quadtree.empty;
                return
            else
                par = qt.parent;
                if strcmp(qt.side, 'nw')
                    node = qt.parent.sw;
                    return
                elseif strcmp(qt.side, 'ne')
                    node = qt.parent.se;
                    return
                end
                mu = south_neighbor(qt.parent);
                if isempty(mu) || isempty(mu.nw)
                    node = mu;
                else
                    if strcmp(qt.side, 'sw')
                        node = mu.nw;
                        return
                    elseif strcmp(qt.side, 'se')
                        node = mu.ne;
                        return
                    end
                end
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function node = east_neighbor(qt)
            if strcmp(qt.side, 'root')
                node = quadtree.empty;
                return
            else
                par = qt.parent;
                if strcmp(qt.side, 'nw')
                    node = qt.parent.ne;
                    return
                elseif strcmp(qt.side, 'sw')
                    node = qt.parent.se;
                    return
                end
                mu = east_neighbor(qt.parent);
                if isempty(mu) || isempty(mu.nw)
                    node = mu;
                else
                    if strcmp(qt.side, 'ne')
                        node = mu.nw;
                        return
                    elseif strcmp(qt.side, 'se')
                        node = mu.sw;
                        return
                    end
                end
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function node = west_neighbor(qt)
            if strcmp(qt.side, 'root')
                node = quadtree.empty;
                return
            else
                par = qt.parent;
                if strcmp(qt.side, 'ne')
                    node = qt.parent.nw;
                    return
                elseif strcmp(qt.side, 'se')
                    node = qt.parent.sw;
                    return
                end
                mu = west_neighbor(qt.parent);
                if isempty(mu) || isempty(mu.nw)
                    node = mu;
                else
                    if strcmp(qt.side, 'nw')
                        node = mu.ne;
                        return
                    elseif strcmp(qt.side, 'sw')
                        node = mu.se;
                        return
                    end
                end
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function find_neighbors(qt)
            qt.nneighbor = qt.north_neighbor();
            if ~isempty(qt.nneighbor) && ~isempty(qt.nneighbor.nw)
                qt.nneighbor = [qt.nneighbor.sw, qt.nneighbor.se];
            end
            
            qt.eneighbor = qt.east_neighbor();
            if ~isempty(qt.eneighbor) && ~isempty(qt.eneighbor.nw)
                qt.eneighbor = [qt.eneighbor.nw, qt.eneighbor.sw];
            end
            
            qt.sneighbor = qt.south_neighbor();
            if ~isempty(qt.sneighbor) && ~isempty(qt.sneighbor.nw)
                qt.sneighbor = [qt.sneighbor.nw, qt.sneighbor.ne];
            end
            
            qt.wneighbor = qt.west_neighbor();
            if ~isempty(qt.wneighbor) && ~isempty(qt.wneighbor.nw)
                qt.wneighbor = [qt.wneighbor.ne, qt.nneighbor.se];
            end
        end
    end
end


