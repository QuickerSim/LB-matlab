classdef quadtreelb < handle
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
        
        unbalanced_leafs(1,:) quadtreelb
        
        xmin
        xmax
        ymin
        ymax
    end
    
    properties
        f
        f_eq
        
        nx
        ny
        
        x
        y
        
        u_ini
        v_ini
        
        u
        v
        
        tau
        omega
        
        w1 = 4/9;
        w2 = 1/9;
        w3 = 1/36;
    end
    
    methods
        function qt = quadtreelb(deep, left, right, bottom, top, mark, parent, side)
            if nargin == 0
                qt.deepness = -1;
                
                qt.xmin = -1;
                qt.xmax = -1;
                qt.ymin = -1;
                qt.ymax = -1;
                
                qt.nw = quadtreelb.empty;
                qt.ne = quadtreelb.empty;
                qt.sw = quadtreelb.empty;
                qt.se = quadtreelb.empty;
                
                qt.marked_for_ref = -1;
                qt.parent = quadtreelb.empty;
                qt.side = "";
                
                qt.nx = -1;
                qt.ny = -1;
                
                qt.tau = -1;
                
                qt.f = [];
            else
                
                qt.deepness = deep;
                
                qt.xmin = left;
                qt.xmax = right;
                qt.ymin = bottom;
                qt.ymax = top;
                
                qt.nw = quadtreelb.empty;
                qt.ne = quadtreelb.empty;
                qt.sw = quadtreelb.empty;
                qt.se = quadtreelb.empty;
                
                qt.marked_for_ref = mark;
                qt.parent = parent;
                qt.side = side;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function create_childrens(qt)
            yhalf = (qt.ymax - qt.ymin)/2;
            xhalf = (qt.xmax - qt.xmin)/2;
            
            qt.nw = quadtreelb(qt.deepness + 1, qt.xmin, qt.xmin + xhalf, qt.ymin + yhalf, qt.ymax, 0, qt, 'nw');
            qt.ne = quadtreelb(qt.deepness + 1, qt.xmin + xhalf, qt.xmax, qt.ymin + yhalf, qt.ymax, 0, qt, 'ne');
            qt.sw = quadtreelb(qt.deepness + 1, qt.xmin, qt.xmin + xhalf, qt.ymin, qt.ymin + yhalf, 0, qt, 'sw');
            qt.se = quadtreelb(qt.deepness + 1, qt.xmin + xhalf, qt.xmax, qt.ymin, qt.ymin + yhalf, 0, qt, 'se');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function find_neighbors(qt)
            qt.nneighbor = qt.north_neighbor();
            if ~isempty(qt.nneighbor) && ~isempty(qt.nneighbor.nw)
                qt.nneighbor = [qt.nneighbor.sw, qt.nneighbor.se];
            end
            
            qt.eneighbor = qt.east_neighbor();
            if ~isempty(qt.eneighbor) && ~isempty(qt.eneighbor.nw)
                qt.eneighbor = [qt.eneighbor.sw, qt.eneighbor.nw];
            end
            
            qt.sneighbor = qt.south_neighbor();
            if ~isempty(qt.sneighbor) && ~isempty(qt.sneighbor.nw)
                qt.sneighbor = [qt.sneighbor.nw, qt.sneighbor.ne];
            end
            
            qt.wneighbor = qt.west_neighbor();
            if ~isempty(qt.wneighbor) && ~isempty(qt.wneighbor.nw)
                qt.wneighbor = [qt.wneighbor.se, qt.wneighbor.ne];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initLB_TestStreamingZeros(qt, nx, ny, tau, uini, vini)
            qt.nx = nx;
            qt.ny = ny;
            
            dx = (qt.xmax - qt.xmin)/nx;
            dy = (qt.ymax - qt.ymin)/ny;
            
            qt.x = linspace(qt.xmin+dx/2, qt.xmax-dx/2, nx);
            qt.y = linspace(qt.ymin+dy/2, qt.ymax-dy/2, ny);
            
%             qt.x = linspace(qt.xmin, qt.xmax, nx);
%             qt.y = linspace(qt.ymin, qt.ymax, ny);
            
            qt.tau = tau;
            qt.omega = 1/tau;
            
            qt.u_ini = uini;
            qt.v_ini = vini;
            
            qt.f = ones(nx+2, ny+2, 9);
            
            qt.f(:,:,1) = 1e-16 * qt.f(:,:,1);
            qt.f(:,:,2) = 1e-16 * qt.f(:,:,2);
            qt.f(:,:,3) = 1e-16 * qt.f(:,:,3);
            qt.f(:,:,4) = 1e-16 * qt.f(:,:,4);
            qt.f(:,:,5) = 1e-16 * qt.f(:,:,5);
            qt.f(:,:,6) = 1e-16 * qt.f(:,:,6);
            qt.f(:,:,7) = 1e-16 * qt.f(:,:,7);
            qt.f(:,:,8) = 1e-16 * qt.f(:,:,8);
            qt.f(:,:,9) = 1e-16 * qt.f(:,:,9);
            
            qt.f_eq = qt.f;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initLB_TestStreaming(qt, nx, ny, tau, uini, vini)
            qt.nx = nx;
            qt.ny = ny;
            
            dx = (qt.xmax - qt.xmin)/nx;
            dy = (qt.ymax - qt.ymin)/ny;
            
            qt.x = linspace(qt.xmin+dx/2, qt.xmax-dx/2, nx);
            qt.y = linspace(qt.ymin+dy/2, qt.ymax-dy/2, ny);
            
            qt.tau = tau;
            qt.omega = 1/tau;
            
            qt.u_ini = uini;
            qt.v_ini = vini;
            
            qt.f = ones(nx+2, ny+2, 9);
            
            qt.f(:,:,1) = qt.w2 * qt.f(:,:,1);
            qt.f(:,:,2) = qt.w3 * qt.f(:,:,2);
            qt.f(:,:,3) = qt.w2 * qt.f(:,:,3);
            qt.f(:,:,4) = qt.w3 * qt.f(:,:,4);
            qt.f(:,:,5) = qt.w2 * qt.f(:,:,5);
            qt.f(:,:,6) = qt.w3 * qt.f(:,:,6);
            qt.f(:,:,7) = qt.w2 * qt.f(:,:,7);
            qt.f(:,:,8) = qt.w3 * qt.f(:,:,8);
            qt.f(:,:,9) = qt.w1 * qt.f(:,:,9);
            
            qt.f(nx/2+1, ny/2+1, 1) = 3;
            qt.f(nx/2+1, ny/2+1, 2) = 3;
            qt.f(nx/2+1, ny/2+1, 3) = 3;
            qt.f(nx/2+1, ny/2+1, 4) = 3;
            qt.f(nx/2+1, ny/2+1, 5) = 3;
            qt.f(nx/2+1, ny/2+1, 6) = 3;
            qt.f(nx/2+1, ny/2+1, 7) = 3;
            qt.f(nx/2+1, ny/2+1, 8) = 3;
            qt.f(nx/2+1, ny/2+1, 9) = 3;
            
            qt.f_eq = qt.f;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initLB(qt, nx, ny, tau, uini, vini)
            qt.nx = nx;
            qt.ny = ny;
            
            dx = (qt.xmax - qt.xmin)/nx;
            dy = (qt.ymax - qt.ymin)/ny;
            
            qt.x = linspace(qt.xmin+dx/2, qt.xmax-dx/2, nx);
            qt.y = linspace(qt.ymin+dy/2, qt.ymax-dy/2, ny);
            
            qt.tau = tau;
            qt.omega = 1/tau;
            
            qt.u_ini = uini;
            qt.v_ini = vini;
            
            qt.f = ones(nx+2, ny+2, 9);
            
            qt.f(:,:,1) = qt.w2 * qt.f(:,:,1);
            qt.f(:,:,2) = qt.w3 * qt.f(:,:,2);
            qt.f(:,:,3) = qt.w2 * qt.f(:,:,3);
            qt.f(:,:,4) = qt.w3 * qt.f(:,:,4);
            qt.f(:,:,5) = qt.w2 * qt.f(:,:,5);
            qt.f(:,:,6) = qt.w3 * qt.f(:,:,6);
            qt.f(:,:,7) = qt.w2 * qt.f(:,:,7);
            qt.f(:,:,8) = qt.w3 * qt.f(:,:,8);
            qt.f(:,:,9) = qt.w1 * qt.f(:,:,9);
            
            qt.f_eq = qt.f;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Collision(qt)
            NX = qt.nx;
            NY = qt.ny;
            
            density = qt.f(2:NX+1, 2:NY+1, 1) + ...
                qt.f(2:NX+1, 2:NY+1, 2) + ...
                qt.f(2:NX+1, 2:NY+1, 3) + ...
                qt.f(2:NX+1, 2:NY+1, 4) + ...
                qt.f(2:NX+1, 2:NY+1, 5) + ...
                qt.f(2:NX+1, 2:NY+1, 6) + ...
                qt.f(2:NX+1, 2:NY+1, 7) + ...
                qt.f(2:NX+1, 2:NY+1, 8) + ...
                qt.f(2:NX+1, 2:NY+1, 9);
            U = (qt.f(2:NX+1, 2:NY+1, 1) + qt.f(2:NX+1, 2:NY+1, 2) + qt.f(2:NX+1, 2:NY+1, 8) - (qt.f(2:NX+1, 2:NY+1, 4) + qt.f(2:NX+1, 2:NY+1, 5) + qt.f(2:NX+1, 2:NY+1, 6)))./density;
            V = (qt.f(2:NX+1, 2:NY+1, 2) + qt.f(2:NX+1, 2:NY+1, 3) + qt.f(2:NX+1, 2:NY+1, 4) - (qt.f(2:NX+1, 2:NY+1, 6) + qt.f(2:NX+1, 2:NY+1, 7) + qt.f(2:NX+1, 2:NY+1, 8)))./density;
            
            qt.u = U;
            qt.v = V;
            
            qt.f_eq(2:NX+1, 2:NY+1, 1) = qt.w2 * density .* (1 + 3*U      + 9/2 * U.^2      - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 3) = qt.w2 * density .* (1 + 3*V      + 9/2 * V.^2      - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 5) = qt.w2 * density .* (1 - 3*U      + 9/2 * U.^2      - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 7) = qt.w2 * density .* (1 - 3*V      + 9/2 * V.^2      - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 2) = qt.w3 * density .* (1 + 3*( U+V) + 9/2 * ( U+V).^2 - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 4) = qt.w3 * density .* (1 + 3*(-U+V) + 9/2 * (-U+V).^2 - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 6) = qt.w3 * density .* (1 - 3*( U+V) + 9/2 * ( U+V).^2 - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 8) = qt.w3 * density .* (1 + 3*( U-V) + 9/2 * ( U-V).^2 - 3/2 * (U.^2 + V.^2));
            qt.f_eq(2:NX+1, 2:NY+1, 9) = qt.w1 * density .* (1                              - 3/2 * (U.^2 + V.^2));
            
            % i=1, Bounceback
            if isempty(qt.wneighbor)
                qt.f(2, 2:NY+1, 1) = qt.f(2, 2:NY+1, 5);
                qt.f(2, 2:NY+1, 2) = qt.f(2, 2:NY+1, 6);
                qt.f(2, 2:NY+1, 8) = qt.f(2, 2:NY+1, 4);
            else
                qt.f(2, 2:NY+1, 1) = qt.omega * qt.f_eq(2, 2:NY+1, 1) + (1-qt.omega) * qt.f(2, 2:NY+1, 1);
                qt.f(2, 2:NY+1, 2) = qt.omega * qt.f_eq(2, 2:NY+1, 2) + (1-qt.omega) * qt.f(2, 2:NY+1, 2);
                qt.f(2, 2:NY+1, 3) = qt.omega * qt.f_eq(2, 2:NY+1, 3) + (1-qt.omega) * qt.f(2, 2:NY+1, 3);
                qt.f(2, 2:NY+1, 4) = qt.omega * qt.f_eq(2, 2:NY+1, 4) + (1-qt.omega) * qt.f(2, 2:NY+1, 4);
                qt.f(2, 2:NY+1, 5) = qt.omega * qt.f_eq(2, 2:NY+1, 5) + (1-qt.omega) * qt.f(2, 2:NY+1, 5);
                qt.f(2, 2:NY+1, 6) = qt.omega * qt.f_eq(2, 2:NY+1, 6) + (1-qt.omega) * qt.f(2, 2:NY+1, 6);
                qt.f(2, 2:NY+1, 7) = qt.omega * qt.f_eq(2, 2:NY+1, 7) + (1-qt.omega) * qt.f(2, 2:NY+1, 7);
                qt.f(2, 2:NY+1, 8) = qt.omega * qt.f_eq(2, 2:NY+1, 8) + (1-qt.omega) * qt.f(2, 2:NY+1, 8);
                qt.f(2, 2:NY+1, 9) = qt.omega * qt.f_eq(2, 2:NY+1, 9) + (1-qt.omega) * qt.f(2, 2:NY+1, 9);
            end
            
            % i=nx, Bounceback
            if isempty(qt.eneighbor)
                qt.f(NX+1, 2:NY+1, 4) = qt.f(NX+1, 2:NY+1, 8);
                qt.f(NX+1, 2:NY+1, 5) = qt.f(NX+1, 2:NY+1, 1);
                qt.f(NX+1, 2:NY+1, 6) = qt.f(NX+1, 2:NY+1, 2);
            else
                qt.f(NX+1, 2:NY+1, 1) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 1) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 1);
                qt.f(NX+1, 2:NY+1, 2) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 2) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 2);
                qt.f(NX+1, 2:NY+1, 3) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 3) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 3);
                qt.f(NX+1, 2:NY+1, 4) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 4) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 4);
                qt.f(NX+1, 2:NY+1, 5) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 5) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 5);
                qt.f(NX+1, 2:NY+1, 6) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 6) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 6);
                qt.f(NX+1, 2:NY+1, 7) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 7) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 7);
                qt.f(NX+1, 2:NY+1, 8) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 8) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 8);
                qt.f(NX+1, 2:NY+1, 9) = qt.omega * qt.f_eq(NX+1, 2:NY+1, 9) + (1-qt.omega) * qt.f(NX+1, 2:NY+1, 9);
            end
            
            % j=1, Bounceback
            if isempty(qt.sneighbor)
                qt.f(2:NX+1, 2, 2) = qt.f(2:NX+1, 2, 6);
                qt.f(2:NX+1, 2, 3) = qt.f(2:NX+1, 2, 7);
                qt.f(2:NX+1, 2, 4) = qt.f(2:NX+1, 2, 8);
            else
                qt.f(2:NX+1, 2, 1) = qt.omega * qt.f_eq(2:NX+1, 2, 1) + (1-qt.omega) * qt.f(2:NX+1, 2, 1);
                qt.f(2:NX+1, 2, 2) = qt.omega * qt.f_eq(2:NX+1, 2, 2) + (1-qt.omega) * qt.f(2:NX+1, 2, 2);
                qt.f(2:NX+1, 2, 3) = qt.omega * qt.f_eq(2:NX+1, 2, 3) + (1-qt.omega) * qt.f(2:NX+1, 2, 3);
                qt.f(2:NX+1, 2, 4) = qt.omega * qt.f_eq(2:NX+1, 2, 4) + (1-qt.omega) * qt.f(2:NX+1, 2, 4);
                qt.f(2:NX+1, 2, 5) = qt.omega * qt.f_eq(2:NX+1, 2, 5) + (1-qt.omega) * qt.f(2:NX+1, 2, 5);
                qt.f(2:NX+1, 2, 6) = qt.omega * qt.f_eq(2:NX+1, 2, 6) + (1-qt.omega) * qt.f(2:NX+1, 2, 6);
                qt.f(2:NX+1, 2, 7) = qt.omega * qt.f_eq(2:NX+1, 2, 7) + (1-qt.omega) * qt.f(2:NX+1, 2, 7);
                qt.f(2:NX+1, 2, 8) = qt.omega * qt.f_eq(2:NX+1, 2, 8) + (1-qt.omega) * qt.f(2:NX+1, 2, 8);
                qt.f(2:NX+1, 2, 9) = qt.omega * qt.f_eq(2:NX+1, 2, 9) + (1-qt.omega) * qt.f(2:NX+1, 2, 9);
            end
            
            % j=NY, Know Velocity
            if isempty(qt.nneighbor)
                densityN = qt.f(2:NX+1, NY+1, 9) + qt.f(2:NX+1, NY+1, 1) + qt.f(2:NX+1, NY+1, 5) + ...
                    2*(qt.f(2:NX+1, NY+1, 3) + qt.f(2:NX+1, NY+1, 4) + qt.f(2:NX+1, NY+1, 2));
                
                qt.f(2:NX+1, NY+1, 6) = qt.f(2:NX+1, NY+1, 2) + 0.5*(qt.f(2:NX+1, NY+1, 1) - qt.f(2:NX+1, NY+1, 5)) - 0.5*qt.u_ini.*densityN + 1/6*densityN*qt.v_ini;
                qt.f(2:NX+1, NY+1, 7) = qt.f(2:NX+1, NY+1, 3) + 2/3*densityN.*qt.v_ini;
                qt.f(2:NX+1, NY+1, 8) = qt.f(2:NX+1, NY+1, 4) + 0.5*(qt.f(2:NX+1, NY+1, 5) - qt.f(2:NX+1, NY+1, 1)) + 0.5*qt.u_ini.*densityN + 1/6*densityN*qt.v_ini;
                
                density = qt.f(2:NX+1, 2:NY+1, 1) + ...
                    qt.f(2:NX+1, 2:NY+1, 2) + ...
                    qt.f(2:NX+1, 2:NY+1, 3) + ...
                    qt.f(2:NX+1, 2:NY+1, 4) + ...
                    qt.f(2:NX+1, 2:NY+1, 5) + ...
                    qt.f(2:NX+1, 2:NY+1, 6) + ...
                    qt.f(2:NX+1, 2:NY+1, 7) + ...
                    qt.f(2:NX+1, 2:NY+1, 8) + ...
                    qt.f(2:NX+1, 2:NY+1, 9);
                U = (qt.f(2:NX+1, 2:NY+1, 1) + qt.f(2:NX+1, 2:NY+1, 2) + qt.f(2:NX+1, 2:NY+1, 8) - (qt.f(2:NX+1, 2:NY+1, 4) + qt.f(2:NX+1, 2:NY+1, 5) + qt.f(2:NX+1, 2:NY+1, 6)))./density;
                V = (qt.f(2:NX+1, 2:NY+1, 2) + qt.f(2:NX+1, 2:NY+1, 3) + qt.f(2:NX+1, 2:NY+1, 4) - (qt.f(2:NX+1, 2:NY+1, 6) + qt.f(2:NX+1, 2:NY+1, 7) + qt.f(2:NX+1, 2:NY+1, 8)))./density;
                
                qt.u = U;
                qt.v = V;
                
                qt.f_eq(2:NX+1, 2:NY+1, 1) = qt.w2 * density .* (1 + 3*U      + 9/2 * U.^2      - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 3) = qt.w2 * density .* (1 + 3*V      + 9/2 * V.^2      - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 5) = qt.w2 * density .* (1 - 3*U      + 9/2 * U.^2      - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 7) = qt.w2 * density .* (1 - 3*V      + 9/2 * V.^2      - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 2) = qt.w3 * density .* (1 + 3*( U+V) + 9/2 * ( U+V).^2 - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 4) = qt.w3 * density .* (1 + 3*(-U+V) + 9/2 * (-U+V).^2 - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 6) = qt.w3 * density .* (1 - 3*( U+V) + 9/2 * ( U+V).^2 - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 8) = qt.w3 * density .* (1 + 3*( U-V) + 9/2 * ( U-V).^2 - 3/2 * (U.^2 + V.^2));
                qt.f_eq(2:NX+1, 2:NY+1, 9) = qt.w1 * density .* (1                              - 3/2 * (U.^2 + V.^2));
                
                qt.f(2:NX+1, NY+1, 1) = qt.omega * qt.f_eq(2:NX+1, NY+1, 1) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 1);
                qt.f(2:NX+1, NY+1, 2) = qt.omega * qt.f_eq(2:NX+1, NY+1, 2) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 2);
                qt.f(2:NX+1, NY+1, 3) = qt.omega * qt.f_eq(2:NX+1, NY+1, 3) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 3);
                qt.f(2:NX+1, NY+1, 4) = qt.omega * qt.f_eq(2:NX+1, NY+1, 4) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 4);
                qt.f(2:NX+1, NY+1, 5) = qt.omega * qt.f_eq(2:NX+1, NY+1, 5) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 5);
                qt.f(2:NX+1, NY+1, 6) = qt.omega * qt.f_eq(2:NX+1, NY+1, 6) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 6);
                qt.f(2:NX+1, NY+1, 7) = qt.omega * qt.f_eq(2:NX+1, NY+1, 7) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 7);
                qt.f(2:NX+1, NY+1, 8) = qt.omega * qt.f_eq(2:NX+1, NY+1, 8) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 8);
                qt.f(2:NX+1, NY+1, 9) = qt.omega * qt.f_eq(2:NX+1, NY+1, 9) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 9);
                
                % j=NY, Bounceback
%                 qt.f(2:NX+1, 2, 2) = qt.f(2:NX+1, 2, 6);
%                 qt.f(2:NX+1, 2, 3) = qt.f(2:NX+1, 2, 7);
%                 qt.f(2:NX+1, 2, 4) = qt.f(2:NX+1, 2, 8);
                
            else
                qt.f(2:NX+1, NY+1, 1) = qt.omega * qt.f_eq(2:NX+1, NY+1, 1) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 1);
                qt.f(2:NX+1, NY+1, 2) = qt.omega * qt.f_eq(2:NX+1, NY+1, 2) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 2);
                qt.f(2:NX+1, NY+1, 3) = qt.omega * qt.f_eq(2:NX+1, NY+1, 3) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 3);
                qt.f(2:NX+1, NY+1, 4) = qt.omega * qt.f_eq(2:NX+1, NY+1, 4) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 4);
                qt.f(2:NX+1, NY+1, 5) = qt.omega * qt.f_eq(2:NX+1, NY+1, 5) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 5);
                qt.f(2:NX+1, NY+1, 6) = qt.omega * qt.f_eq(2:NX+1, NY+1, 6) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 6);
                qt.f(2:NX+1, NY+1, 7) = qt.omega * qt.f_eq(2:NX+1, NY+1, 7) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 7);
                qt.f(2:NX+1, NY+1, 8) = qt.omega * qt.f_eq(2:NX+1, NY+1, 8) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 8);
                qt.f(2:NX+1, NY+1, 9) = qt.omega * qt.f_eq(2:NX+1, NY+1, 9) + (1-qt.omega) * qt.f(2:NX+1, NY+1, 9);
            end
            
            qt.f(3:NX, 3:NY, 1) = qt.omega * qt.f_eq(3:NX, 3:NY, 1) + (1-qt.omega) * qt.f(3:NX, 3:NY, 1);
            qt.f(3:NX, 3:NY, 2) = qt.omega * qt.f_eq(3:NX, 3:NY, 2) + (1-qt.omega) * qt.f(3:NX, 3:NY, 2);
            qt.f(3:NX, 3:NY, 3) = qt.omega * qt.f_eq(3:NX, 3:NY, 3) + (1-qt.omega) * qt.f(3:NX, 3:NY, 3);
            qt.f(3:NX, 3:NY, 4) = qt.omega * qt.f_eq(3:NX, 3:NY, 4) + (1-qt.omega) * qt.f(3:NX, 3:NY, 4);
            qt.f(3:NX, 3:NY, 5) = qt.omega * qt.f_eq(3:NX, 3:NY, 5) + (1-qt.omega) * qt.f(3:NX, 3:NY, 5);
            qt.f(3:NX, 3:NY, 6) = qt.omega * qt.f_eq(3:NX, 3:NY, 6) + (1-qt.omega) * qt.f(3:NX, 3:NY, 6);
            qt.f(3:NX, 3:NY, 7) = qt.omega * qt.f_eq(3:NX, 3:NY, 7) + (1-qt.omega) * qt.f(3:NX, 3:NY, 7);
            qt.f(3:NX, 3:NY, 8) = qt.omega * qt.f_eq(3:NX, 3:NY, 8) + (1-qt.omega) * qt.f(3:NX, 3:NY, 8);
            qt.f(3:NX, 3:NY, 9) = qt.omega * qt.f_eq(3:NX, 3:NY, 9) + (1-qt.omega) * qt.f(3:NX, 3:NY, 9);
            
%             qt.f([2 NX+1], [2 NY+1], 1) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 1) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 1);
%             qt.f([2 NX+1], [2 NY+1], 2) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 2) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 2);
%             qt.f([2 NX+1], [2 NY+1], 3) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 3) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 3);
%             qt.f([2 NX+1], [2 NY+1], 4) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 4) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 4);
%             qt.f([2 NX+1], [2 NY+1], 5) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 5) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 5);
%             qt.f([2 NX+1], [2 NY+1], 6) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 6) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 6);
%             qt.f([2 NX+1], [2 NY+1], 7) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 7) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 7);
%             qt.f([2 NX+1], [2 NY+1], 8) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 8) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 8);
%             qt.f([2 NX+1], [2 NY+1], 9) = qt.omega * qt.f_eq([2 NX+1], [2 NY+1], 9) + (1-qt.omega) * qt.f([2 NX+1], [2 NY+1], 9);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Streaming(qt)
            NX = qt.nx;
            NY = qt.ny;
            
            % interior
            qt.f(2:NX+1, 2:NY+1, 1) = qt.f(1:NX  , 2:NY+1, 1);
            qt.f(2:NX+1, 2:NY+1, 2) = qt.f(1:NX  , 1:NY  , 2);
            qt.f(2:NX+1, 2:NY+1, 3) = qt.f(2:NX+1, 1:NY  , 3);
            qt.f(2:NX+1, 2:NY+1, 4) = qt.f(3:NX+2, 1:NY  , 4);
            qt.f(2:NX+1, 2:NY+1, 5) = qt.f(3:NX+2, 2:NY+1, 5);
            qt.f(2:NX+1, 2:NY+1, 6) = qt.f(3:NX+2, 3:NY+2, 6);
            qt.f(2:NX+1, 2:NY+1, 7) = qt.f(2:NX+1, 3:NY+2, 7);
            qt.f(2:NX+1, 2:NY+1, 8) = qt.f(1:NX  , 3:NY+2, 8);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GhostInfo(qt)
            NX = qt.nx;
            NY = qt.ny;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% EAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.eneighbor) == 1 && qt.deepness == qt.eneighbor.deepness
                qt.f(NX+2, 2:NY+1, 5) = qt.eneighbor.f(2, 2:NY+1, 5);
                qt.f(NX+2, 2:NY+1, 4) = qt.eneighbor.f(2, 2:NY+1, 4);
                qt.f(NX+2, 2:NY+1, 6) = qt.eneighbor.f(2, 2:NY+1, 6);
            end
            
%             if length(qt.eneighbor.nneighbor) == 1
%                 qt.f(NX+2, NY+2, 6) = qt.eneighbor.nneighbor.f(2, 2, 6);
%             end
% 
%             if length(qt.eneighbor.sneighbor) == 1
%                 qt.f(NX+2, 1, 4) = qt.eneighbor.sneighbor.f(2, NY+1, 4);
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% NORTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.nneighbor) == 1 && qt.deepness == qt.nneighbor.deepness
                qt.f(2:NX+1, NY+2, 7) = qt.nneighbor.f(2:NX+1, 2, 7);
                qt.f(2:NX+1, NY+2, 6) = qt.nneighbor.f(2:NX+1, 2, 6);
                qt.f(2:NX+1, NY+2, 8) = qt.nneighbor.f(2:NX+1, 2, 8);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% WEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.wneighbor) == 1 && qt.deepness == qt.wneighbor.deepness
                qt.f(1, 2:NY+1, 1) = qt.wneighbor.f(NX+1, 2:NY+1, 1);
                qt.f(1, 2:NY+1, 2) = qt.wneighbor.f(NX+1, 2:NY+1, 2);
                qt.f(1, 2:NY+1, 8) = qt.wneighbor.f(NX+1, 2:NY+1, 8);
            end
            
%             if length(qt.wneighbor.nneighbor) == 1
%                 qt.f(1, NY+2, 8) = qt.wneighbor.nneighbor.f(NX+1, 2, 8);
%             end
% 
%             if length(qt.wneighbor.sneighbor) == 1
%                 qt.f(1, 1, 2) = qt.wneighbor.sneighbor.f(NX+1, NY+1, 2);
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% SOUTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.sneighbor) == 1 && qt.deepness == qt.sneighbor.deepness    
                qt.f(2:NX+1, 1, 3) = qt.sneighbor.f(2:NX+1, NY+1, 3);
                qt.f(2:NX+1, 1, 2) = qt.sneighbor.f(2:NX+1, NY+1, 2);
                qt.f(2:NX+1, 1, 4) = qt.sneighbor.f(2:NX+1, NY+1, 4);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InterpolateDown(qt)
            NX = qt.nx;
            NY = qt.ny;
            
            %%%%%%%%%%%%%%%%%%%% EAST %%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.eneighbor) == 1 && qt.deepness > qt.eneighbor.deepness
                if strcmp(qt.side, 'se')
                    qt.f(NX+2, 2:NY+1, 5) = qt.eneighbor.f(2, repelem(2:NY/2+1,2), 5);
                    qt.f(NX+2, 2:NY+1, 4) = qt.eneighbor.f(2, repelem(2:NY/2+1,2), 4);
                    qt.f(NX+2, 2:NY+1, 6) = qt.eneighbor.f(2, repelem(2:NY/2+1,2), 6);
                elseif strcmp(qt.side, 'ne')
                    qt.f(NX+2, 2:NY+1, 5) = qt.eneighbor.f(2, repelem((NY+2)/2+1:NY+1,2), 5);
                    qt.f(NX+2, 2:NY+1, 4) = qt.eneighbor.f(2, repelem((NY+2)/2+1:NY+1,2), 4);
                    qt.f(NX+2, 2:NY+1, 6) = qt.eneighbor.f(2, repelem((NY+2)/2+1:NY+1,2), 6); 
                end
            end
            %%%%%%%%%%%%%%%%%%%% NORTH %%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.nneighbor) == 1 && qt.deepness > qt.nneighbor.deepness
                if strcmp(qt.side, 'nw')
                    qt.f(2:NX+1, NY+2, 7) = qt.nneighbor.f(repelem(2:NX/2+1, 2), 2, 7);
                    qt.f(2:NX+1, NY+2, 6) = qt.nneighbor.f(repelem(2:NX/2+1, 2), 2, 6);
                    qt.f(2:NX+1, NY+2, 8) = qt.nneighbor.f(repelem(2:NX/2+1, 2), 2, 8);
                elseif strcmp(qt.side, 'ne')
                    qt.f(2:NX+1, NY+2, 7) = qt.nneighbor.f(repelem((NX+2)/2+1:NX+1, 2), 2, 7);
                    qt.f(2:NX+1, NY+2, 6) = qt.nneighbor.f(repelem((NX+2)/2+1:NX+1, 2), 2, 6);
                    qt.f(2:NX+1, NY+2, 8) = qt.nneighbor.f(repelem((NX+2)/2+1:NX+1, 2), 2, 8);
                end
            end
            %%%%%%%%%%%%%%%%%%%%% WEST %%%%%%%%%%%%%%%%%%%%%%
            if length(qt.wneighbor) == 1 && qt.deepness > qt.wneighbor.deepness
                if strcmp(qt.side, 'sw')
                    qt.f(1, 2:NY+1, 1) = qt.wneighbor.f(NX+1, repelem(2:NY/2+1,2), 1);
                    qt.f(1, 2:NY+1, 2) = qt.wneighbor.f(NX+1, repelem(2:NY/2+1,2), 2);
                    qt.f(1, 2:NY+1, 8) = qt.wneighbor.f(NX+1, repelem(2:NY/2+1,2), 8);
                elseif strcmp(qt.side, 'nw')
                    qt.f(1, 2:NY+1, 1) = qt.wneighbor.f(NX+1, repelem((NY+2)/2+1:NY+1,2), 1);
                    qt.f(1, 2:NY+1, 2) = qt.wneighbor.f(NX+1, repelem((NY+2)/2+1:NY+1,2), 2);
                    qt.f(1, 2:NY+1, 8) = qt.wneighbor.f(NX+1, repelem((NY+2)/2+1:NY+1,2), 8);
                end
            end
            %%%%%%%%%%%%%%%%%%%%% SOUTH %%%%%%%%%%%%%%%%%%%%%%
            if length(qt.sneighbor) == 1 && qt.deepness > qt.sneighbor.deepness
                if strcmp(qt.side, 'sw')
                    qt.f(2:NX+1, 1, 3) = qt.sneighbor.f(repelem(2:NX/2+1, 2), NY+1, 3);
                    qt.f(2:NX+1, 1, 2) = qt.sneighbor.f(repelem(2:NX/2+1, 2), NY+1, 2);
                    qt.f(2:NX+1, 1, 4) = qt.sneighbor.f(repelem(2:NX/2+1, 2), NY+1, 4);
                elseif strcmp(qt.side, 'se')
                    qt.f(2:NX+1, 1, 3) = qt.sneighbor.f(repelem((NX+2)/2+1:NX+1, 2), NY+1, 3);
                    qt.f(2:NX+1, 1, 2) = qt.sneighbor.f(repelem((NX+2)/2+1:NX+1, 2), NY+1, 2);
                    qt.f(2:NX+1, 1, 4) = qt.sneighbor.f(repelem((NX+2)/2+1:NX+1, 2), NY+1, 4);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InterpolateUp(qt)
            NX = qt.nx;
            NY = qt.ny;
            
            %%%%%%%%%%%%%%%%%%% SOUTH %%%%%%%%%%%%%%%%%%%%%%
            if length(qt.sneighbor) == 2
                qt.f(2:NX+1, 1, 3) = [(qt.sneighbor(1).f(2:2:NX+1, NY+1, 3)+qt.sneighbor(1).f(3:2:NX+1, NY+1, 3))/2; ...
                                      (qt.sneighbor(2).f(2:2:NX+1, NY+1, 3)+qt.sneighbor(2).f(3:2:NX+1, NY+1, 3))/2];

                qt.f(2:NX+1, 1, 2) = [(qt.sneighbor(1).f(2:2:NX+1, NY+1, 2)+qt.sneighbor(1).f(3:2:NX+1, NY+1, 2))/2; ...
                                      (qt.sneighbor(2).f(2:2:NX+1, NY+1, 2)+qt.sneighbor(2).f(3:2:NX+1, NY+1, 2))/2];

                qt.f(2:NX+1, 1, 4) = [(qt.sneighbor(1).f(2:2:NX+1, NY+1, 4)+qt.sneighbor(1).f(3:2:NX+1, NY+1, 4))/2; ...
                                      (qt.sneighbor(2).f(2:2:NX+1, NY+1, 4)+qt.sneighbor(2).f(3:2:NX+1, NY+1, 4))/2];
            end
            %%%%%%%%%%%%%%%%%%%%% WEST %%%%%%%%%%%%%%%%%%%%%%
            if length(qt.wneighbor) == 2
                qt.f(1, 2:NY+1, 1) = [(qt.wneighbor(1).f(NX+1, 2:2:NY+1, 1)+qt.wneighbor(1).f(NX+1, 3:2:NY+1, 1))/2 ...
                                      (qt.wneighbor(2).f(NX+1, 2:2:NY+1, 1)+qt.wneighbor(2).f(NX+1, 3:2:NY+1, 1))/2];

                qt.f(1, 2:NY+1, 2) = [(qt.wneighbor(1).f(NX+1, 2:2:NY+1, 2)+qt.wneighbor(1).f(NX+1, 3:2:NY+1, 2))/2 ...
                                      (qt.wneighbor(2).f(NX+1, 2:2:NY+1, 2)+qt.wneighbor(2).f(NX+1, 3:2:NY+1, 2))/2];

                qt.f(1, 2:NY+1, 8) = [(qt.wneighbor(1).f(NX+1, 2:2:NY+1, 8)+qt.wneighbor(1).f(NX+1, 3:2:NY+1, 8))/2 ...
                                      (qt.wneighbor(2).f(NX+1, 2:2:NY+1, 8)+qt.wneighbor(2).f(NX+1, 3:2:NY+1, 8))/2];
            end
            %%%%%%%%%%%%%%%%%%%% NORTH %%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.nneighbor) == 2
                qt.f(2:NX+1, NY+2, 7) = [(qt.nneighbor(1).f(2:2:NX+1, 2, 7)+qt.nneighbor(1).f(3:2:NX+1, 2, 7))/2; ...
                                         (qt.nneighbor(2).f(2:2:NX+1, 2, 7)+qt.nneighbor(2).f(3:2:NX+1, 2, 7))/2];

                qt.f(2:NX+1, NY+2, 6) = [(qt.nneighbor(1).f(2:2:NX+1, 2, 6)+qt.nneighbor(1).f(3:2:NX+1, 2, 6))/2; ...
                                         (qt.nneighbor(2).f(2:2:NX+1, 2, 6)+qt.nneighbor(2).f(3:2:NX+1, 2, 6))/2];

                qt.f(2:NX+1, NY+2, 8) = [(qt.nneighbor(1).f(2:2:NX+1, 2, 8)+qt.nneighbor(1).f(3:2:NX+1, 2, 8))/2; ...
                                         (qt.nneighbor(2).f(2:2:NX+1, 2, 8)+qt.nneighbor(2).f(3:2:NX+1, 2, 8))/2];
            end
            %%%%%%%%%%%%%%%%%%%% EAST %%%%%%%%%%%%%%%%%%%%%%%%%
            if length(qt.eneighbor) == 2
                qt.f(NX+2, 2:NY+1, 5) = [(qt.eneighbor(1).f(2, 2:2:NY+1, 5)+qt.eneighbor(1).f(2, 3:2:NY+1, 5))/2 ...
                                         (qt.eneighbor(2).f(2, 2:2:NY+1, 5)+qt.eneighbor(2).f(2, 3:2:NY+1, 5))/2];

                qt.f(NX+2, 2:NY+1, 4) = [(qt.eneighbor(1).f(2, 2:2:NY+1, 4)+qt.eneighbor(1).f(2, 3:2:NY+1, 4))/2 ...
                                         (qt.eneighbor(2).f(2, 2:2:NY+1, 4)+qt.eneighbor(2).f(2, 3:2:NY+1, 4))/2];

                qt.f(NX+2, 2:NY+1, 6) = [(qt.eneighbor(1).f(2, 2:2:NY+1, 6)+qt.eneighbor(1).f(2, 3:2:NY+1, 6))/2 ...
                                         (qt.eneighbor(2).f(2, 2:2:NY+1, 6)+qt.eneighbor(2).f(2, 3:2:NY+1, 6))/2];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function drawLB(qt)
            imagesc(qt.x, qt.y, sqrt(qt.u'.^2 + qt.v'.^2))
            %contour(qt.x, qt.y, sqrt(qt.u'.^2 + qt.v'.^2), [0:0.001:0.5]);
            colormap(jet(6400))
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end















