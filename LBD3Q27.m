classdef LBD3Q27 < handle
    % Opis oznaczeń
    % X - x coordinate in the positive direction
    % x - x coordinate in the negative direction
    % Y - y coordinate in the positive direction
    % y - y coordinate in the negative direction
    % Z - z coordinate in the positive direction
    % z - z coordinate in the negative direction
    
    % Opis kierunków
    %  1 - X
    %  2 - Y
    %  3 - Z
    %  4 - x
    %  5 - y
    %  6 - z
    %  7 - YZ
    %  8 - Yz
    %  9 - yZ
    % 10 - yz
    % 11 - XZ
    % 12 - xZ
    % 13 - Xz
    % 14 - xz
    % 15 - XY
    % 16 - Xy
    % 17 - xY
    % 18 - xy
    % 19 - XYZ
    % 20 - xYZ
    % 21 - XyZ
    % 22 - xyZ
    % 23 - XYz
    % 24 - xYz
    % 25 - Xyz
    % 26 - xyz
    % 27 - center
    properties (Access = private)
        w1  = 2/27;
        w2  = 2/27;
        w3  = 2/27;
        w4  = 2/27;
        w5  = 2/27;
        w6  = 2/27;
        
        w7  = 1/54;
        w8  = 1/54;
        w9  = 1/54;
        w10 = 1/54;
        w11 = 1/54;
        w12 = 1/54;
        w13 = 1/54;
        w14 = 1/54;
        w15 = 1/54;
        w16 = 1/54;
        w17 = 1/54;
        w18 = 1/54;
        
        w19 = 1/126;
        w20 = 1/126;
        w21 = 1/126;
        w22 = 1/126;
        w23 = 1/126;
        w24 = 1/126;
        w25 = 1/126;
        w26 = 1/126;
        
        w27 = 8/27;
    end
    
    properties
        N % 3x1 vector - size of subgrid in every direction
        k % suize of ghosts node thickness
        % (mostly n(1)=n(2)=n(3))
        f % (N+k*2)x(N+k*2)x(N+k*2)x27 - populations in 27 directions for
        % 3D where 'k' is size of ghost layer- can be generated by template or script
        treeID % ID of node (leaf) in octree (global numeration)
        layer % number of layer on which node (leaf) exist
        coords % can be computed, but for better performence its saved
        vel
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = LBD3Q27(treeID, N, k)
            if nargin == 3
                self.treeID = treeID;
                
                if treeID == 0
                    L = 0;
                elseif (treeID >= 1 && treeID <= 8)
                    L = 1;
                elseif (treeID >= 9 && treeID <= 72)
                    L = 2;
                elseif (treeID >= 73 && treeID <= 584)
                    L = 3;
                elseif (treeID >= 585 && treeID <= 4680)
                    L = 4;
                elseif (treeID >= 4681 && treeID <= 37448)
                    L = 5;
                elseif (treeID >= 37449 && treeID <= 299592)
                    L = 6;
                elseif (treeID >= 299593 && treeID <= 2396744)
                    L = 7;
                elseif (treeID >= 2396745 && treeID <= 19173960)
                    L = 8;
                elseif (treeID >= 19173961 && treeID <= 153391688)
                    L = 9;
                elseif (treeID >= 153391689 && treeID <= 1227133512)
                    L = 10;
                elseif (treeID >= 1227133513 && treeID <= 9817068104)
                    L = 11;
                elseif (treeID >= 9817068105 && treeID <= 78536544840)
                    L = 12;
                elseif (treeID >= 78536544841 && treeID <= 628292358728)
                    L = 13;
                elseif (treeID >= 628292358729 && treeID <= 5026338869833)
                    L = 14;
                else
                    error('Too deep layer, not implemented')
                end
                
                self.layer = L;
                self.N = [N N N];
                self.k = k;
                self.f = ones(N+2*k, N+2*k, N+2*k, 27);
                
                self.f(:,:,:, 1) = self.w1;
                self.f(:,:,:, 2) = self.w2;
                self.f(:,:,:, 3) = self.w3;
                self.f(:,:,:, 4) = self.w4;
                self.f(:,:,:, 5) = self.w5;
                self.f(:,:,:, 6) = self.w6;
                self.f(:,:,:, 7) = self.w7;
                self.f(:,:,:, 8) = self.w8;
                self.f(:,:,:, 9) = self.w9;
                self.f(:,:,:,10) = self.w10;
                self.f(:,:,:,11) = self.w11;
                self.f(:,:,:,12) = self.w12;
                self.f(:,:,:,13) = self.w13;
                self.f(:,:,:,14) = self.w14;
                self.f(:,:,:,15) = self.w15;
                self.f(:,:,:,16) = self.w16;
                self.f(:,:,:,17) = self.w17;
                self.f(:,:,:,18) = self.w18;
                self.f(:,:,:,19) = self.w19;
                self.f(:,:,:,20) = self.w20;
                self.f(:,:,:,21) = self.w21;
                self.f(:,:,:,22) = self.w22;
                self.f(:,:,:,23) = self.w23;
                self.f(:,:,:,24) = self.w24;
                self.f(:,:,:,25) = self.w25;
                self.f(:,:,:,26) = self.w26;
                
%                 self.coords = self.GlobalToCoords();
            elseif nargin == 0
                self.treeID = [];
                self.layer = [];
                self.N = [];
                self.k = [];
                self.f = [];
                self.coords = [];
            else
                error('Not enough input arguments. KJ');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function coords = GlobalToCoords(self)
            dim = 3;
            
            elemnums = cumsum((2^dim).^[0:15]);
            
            for k = 3:-1:1
                coords(k) = bin2dec(join(string(fliplr(bitget(self.treeID - ...
                elemnums(self.layer), k:3:25))),''));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Collision(self, omega)
            NX = self.N(1);
            NY = self.N(2);
            NZ = self.N(3);
            ki = self.k;
            
            
            density = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 1) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 2) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 3) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 4) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 5) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 6) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 7) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 8) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 9) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) + ...
                      self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 27);
            
            u = (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  1) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) - ...
                (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  4) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) ...
                ) ...
                )./density;
            
            v = (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  2) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  7) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  8) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) - ...
                (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  5) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  9) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) ...
                ) ...
                )./density;
            
            w = (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  3) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  7) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  9) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) - ...
                (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  6) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  8) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) ...
                ) ...
                )./density;
            
            self.vel{1} = u;
            self.vel{2} = v;
            self.vel{3} = w;
            
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  1) = omega * self.w1 * density .* (1 + 3 * u + 9/2 * u.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 1);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  2) = omega * self.w2 * density .* (1 + 3 * v + 9/2 * v.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 2);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  3) = omega * self.w3 * density .* (1 + 3 * w + 9/2 * w.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 3);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  4) = omega * self.w4 * density .* (1 - 3 * u + 9/2 * u.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 4);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  5) = omega * self.w5 * density .* (1 - 3 * v + 9/2 * v.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 5);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  6) = omega * self.w6 * density .* (1 - 3 * w + 9/2 * w.^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 6);
            
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  7) = omega * self.w7 * density .* (1 + 3 * ( v+w) + 9/2 * ( v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 7);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  8) = omega * self.w8 * density .* (1 + 3 * ( v-w) + 9/2 * ( v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 8);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  9) = omega * self.w9 * density .* (1 + 3 * (-v+w) + 9/2 * (-v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 9);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10) = omega * self.w10 * density .* (1 + 3 * (-v-w) + 9/2 * (-v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) = omega * self.w11 * density .* (1 + 3 * ( u+w) + 9/2 * ( u+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) = omega * self.w12 * density .* (1 + 3 * (-u+w) + 9/2 * (-u+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) = omega * self.w13 * density .* (1 + 3 * ( u-w) + 9/2 * ( u-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) = omega * self.w14 * density .* (1 + 3 * (-u-w) + 9/2 * (-u-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) = omega * self.w15 * density .* (1 + 3 * ( u+v) + 9/2 * ( u+v).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) = omega * self.w16 * density .* (1 + 3 * ( u-v) + 9/2 * ( u-v).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17) = omega * self.w17 * density .* (1 + 3 * (-u+v) + 9/2 * (-u+v).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) = omega * self.w18 * density .* (1 + 3 * (-u-v) + 9/2 * (-u-v).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18);
            
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) = omega * self.w19 * density .* (1 + 3 * ( u+v+w) + 9/2 * ( u+v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) = omega * self.w20 * density .* (1 + 3 * (-u+v+w) + 9/2 * (-u+v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) = omega * self.w21 * density .* (1 + 3 * ( u-v+w) + 9/2 * ( u-v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) = omega * self.w22 * density .* (1 + 3 * (-u-v+w) + 9/2 * (-u-v+w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) = omega * self.w23 * density .* (1 + 3 * ( u+v-w) + 9/2 * ( u+v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) = omega * self.w24 * density .* (1 + 3 * (-u+v-w) + 9/2 * (-u+v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) = omega * self.w25 * density .* (1 + 3 * ( u-v-w) + 9/2 * ( u-v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) = omega * self.w26 * density .* (1 + 3 * (-u-v-w) + 9/2 * (-u-v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26);
            
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 27) = omega * self.w27 * density .* (1 + 3 * (-u-v-w) + 9/2 * (-u-v-w).^2 - 3/2 * (u.^2 + v.^2 + w.^2) ) + ...
                (1-omega) * self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 27);
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Streaming(self)
            NX = self.N(1);
            NY = self.N(2);
            NZ = self.N(3);
            ki = self.k;
            % interior
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  1) = self.f((ki):(NX+ki-1)  , (ki+1):(NY+ki)  , (ki+1):(NZ+ki)  , 1);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  2) = self.f((ki+1):(NX+ki)  , (ki):(NY+ki-1)  , (ki+1):(NZ+ki)  , 2);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  3) = self.f((ki+1):(NX+ki)  , (ki+1):(NY+ki)  , (ki):(NZ+ki-1)  , 3);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  4) = self.f((ki+2):(NX+ki+1), (ki+1):(NY+ki)  , (ki+1):(NZ+ki)  , 4);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  5) = self.f((ki+1):(NX+ki)  , (ki+2):(NY+ki+1), (ki+1):(NZ+ki)  , 5);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  6) = self.f((ki+1):(NX+ki)  , (ki+1):(NY+ki)  , (ki+2):(NZ+ki+1), 6);

            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  7) = self.f((ki+1):(NX+ki)  , (ki):(NY+ki-1)  , (ki):(NZ+ki-1)  , 7);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  8) = self.f((ki+1):(NX+ki)  , (ki):(NY+ki-1)  , (ki+2):(NZ+ki+1), 8);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki),  9) = self.f((ki+1):(NX+ki)  , (ki+2):(NY+ki+1), (ki):(NZ+ki-1)  , 9);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 10) = self.f((ki+1):(NX+ki)  , (ki+2):(NY+ki+1), (ki+2):(NZ+ki+1), 10);

            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) = self.f((ki):(NX+ki-1)  , (ki+1):(NY+ki)  , (ki):(NZ+ki-1)  , 11);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) = self.f((ki+2):(NX+ki+1), (ki+1):(NY+ki)  , (ki):(NZ+ki-1)  , 12);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) = self.f((ki):(NX+ki-1)  , (ki+1):(NY+ki)  , (ki+2):(NZ+ki+1), 13);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) = self.f((ki+2):(NX+ki+1), (ki+1):(NY+ki)  , (ki+2):(NZ+ki+1), 14);

            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) = self.f((ki):(NX+ki-1)  , (ki):(NY+ki-1)  , (ki+1):(NZ+ki)  , 15);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) = self.f((ki):(NX+ki-1)  , (ki+2):(NY+ki+1), (ki+1):(NZ+ki)  , 16);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 17) = self.f((ki+2):(NX+ki+1), (ki):(NY+ki-1)  , (ki+1):(NZ+ki)  , 17);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) = self.f((ki+2):(NX+ki+1), (ki+2):(NY+ki+1), (ki+1):(NZ+ki)  , 18);

            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) = self.f((ki):(NX+ki-1)  , (ki):(NY+ki-1)  , (ki):(NZ+ki-1)  , 19);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) = self.f((ki+2):(NX+ki+1), (ki):(NY+ki-1)  , (ki):(NZ+ki-1)  , 20);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) = self.f((ki):(NX+ki-1)  , (ki+2):(NY+ki+1), (ki):(NZ+ki-1)  , 21);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) = self.f((ki+2):(NX+ki+1), (ki+2):(NY+ki+1), (ki):(NZ+ki-1)  , 22);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) = self.f((ki):(NX+ki-1)  , (ki):(NY+ki-1)  , (ki+2):(NZ+ki+1), 23);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) = self.f((ki+2):(NX+ki+1), (ki):(NY+ki-1)  , (ki+2):(NZ+ki+1), 24);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) = self.f((ki):(NX+ki-1)  , (ki+2):(NY+ki+1)  , (ki+2):(NZ+ki+1), 25);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) = self.f((ki+2):(NX+ki+1), (ki+2):(NY+ki+1)  , (ki+2):(NZ+ki+1), 26);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function BC(self)
            ki = self.k;
            % x
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 21) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 24);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 11) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 14);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 19) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 26);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 27);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki),  1) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki),  4);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 15) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 18);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 25) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 20);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 13) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 12);
            self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 23) = self.f(ki+1, (ki+1):(NY+ki), (ki+1):(NZ+ki), 22);
            
            % X
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 24) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 21)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 14) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 11)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 26) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 19)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 27) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 16) 
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki),  4) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki),  1)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 18) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 15)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 20) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 25)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 12) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 13)
            self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 22) = self.f(NX+ki, (ki+1):(NY+ki), (ki+1):(NZ+ki), 23)
            
            % y
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 19) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 26);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki),  7) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 10);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 20) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 25);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 15) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 18);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki),  2) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki),  5);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 17) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 16);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 23) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 22);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki),  8) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki),  9);
            self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 24) = self.f((ki+1):(NX+ki), ki+1, (ki+1):(NZ+ki), 21);
            
            % Y
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 26) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 19);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 10) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki),  7);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 25) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 20);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 18) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 15);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki),  5) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki),  2);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 16) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 17);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 22) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 23);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki),  9) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki),  8);
            self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 21) = self.f((ki+1):(NX+ki), NY+ki, (ki+1):(NZ+ki), 24);
            
            % z
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 21) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 24);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1,  9) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1,  8);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 22) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 23);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 11) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 14);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1,  3) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1,  6);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 12) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 13);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 19) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 26);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1,  7) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 10);
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 20) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), ki+1, 25);
            
            % Z - top
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 24) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 21);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  8) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  9);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 23) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 22);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 14) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 11);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  6) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  3);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 13) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 12);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 26) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 19);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 10) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  7);
            % self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 25) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 20);

            densityN = ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 16) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  1) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 15) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  5) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 27) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  2) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 18) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  4) + ...
                self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 17) + ...
               2*(self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 21) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 11) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 19) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  9) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  3) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  7) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 22) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 12) + ...
                  self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 20) ...
                 );

            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 25) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 2) + ...
                0.5 * (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 1) - self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 5)) - 0.5*qt.u_ini.*densityN + 1/6*densityN*qt.v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 13) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 23) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 10) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  6) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki,  8) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 26) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 14) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 3) + ... 
                2/3 * densityN .* v_ini;
            self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 24) = self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 4) + ... 
                0.5 * (self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 5) - self.f((ki+1):(NX+ki), (ki+1):(NY+ki), NZ+ki, 1)) + 0.5*qt.u_ini.*densityN + 1/6*densityN*qt.v_ini;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end































