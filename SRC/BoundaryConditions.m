function [FCond, BC_fixed, BC] = BoundaryConditions(nelx, nely, XYZ, LE, ndof, BVP, ForceMagnitude, theta_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BVP ==0 %% uniaxial
    BC_fixed = 1:(nely+1)*ndof;
    % point of interest to be investigated
    top_right_node = (nelx+1)*(nely+1); % nodes at which force is applied
    bottom_right_node = (nelx)*(nely+1)+1;
    Fnodes = [bottom_right_node:top_right_node];
    Fdir = 1*ones(size(Fnodes)); % direction of force
    Ftheta = theta_0*ones(size(Fnodes)); % theta of the applied force
    F0 = ForceMagnitude; % magnitude of force
    Fvals = F0*ones(size(Fnodes));
    POI = bottom_right_node+2;
end


if BVP ==1 %% cantilever beam at top right
    BC_fixed = 1:(nely+1)*ndof;
    % point of interest to be investigated
    poi = (nelx+1)*(nely+1); % at which force is applied
    Fnodes = poi;
    Fdir = -2*ones(size(Fnodes)); % direction of force
    Ftheta = theta_0*ones(size(Fnodes));
    F0 = ForceMagnitude; % magnitude of force
    Fvals = F0*ones(size(Fnodes));
end

if BVP ==101 %% cantilever beam at top right
    BC_fixed = 1:(nely+1)*ndof;
    % point of interest to be investigated
    poi = (nelx+1)*(nely+1) - (nely+1) +1; % at which force is applied
    Fnodes = poi;
    Fdir = -2*ones(size(Fnodes)); % direction of force
    Ftheta = theta_0*ones(size(Fnodes));
    F0 = ForceMagnitude; % magnitude of force
    Fvals = F0*ones(size(Fnodes));
end

%% MBB beam
if BVP ==2

    node_lower_corner_right = ((nelx+1)*(nely+1));

    node_lower_corner_left = nely+1;

    BC_fixed = [node_lower_corner_left*2-1 node_lower_corner_left*2  ...
        node_lower_corner_right*2-1 node_lower_corner_right*2];

    % force nodes
    if mod((nelx+1),2) == 0
        Fnodes = [(nely+1)*ceil((nelx+1)/2) (nely+1)*ceil((nelx+1)/2)+nely+1];
    end
    if mod((nelx+1),2) ~=0
        Fnodes = (nely+1)*ceil((nelx+1)/2);
    end
    poi = Fnodes;
    Fdir = -2*ones(size(Fnodes)); % direction of force
    Ftheta = theta_0*ones(size(Fnodes));

    F0 = ForceMagnitude; % magnitude of force
    Fvals = F0*ones(size(Fnodes));

end

if BVP == 4 %% cantilever beam at middle right

    if rem (nely, 2)~=0
        error("nely must be even")
    end
    BC_fixed = 1:(nely+1)*ndof;
    % point of interest to be investigated
    poi = (nely+1)*(nelx) + 0.5*(nely) +1; % at which node force is applied
    Fnodes = poi;
    Fdir = -2*ones(size(Fnodes)); % direction of force
    Ftheta = theta_0*ones(size(Fnodes));
    F0 = ForceMagnitude; % magnitude of force
    Fvals = F0*ones(size(Fnodes));
end



FCond = [Fnodes' Fdir' Fvals' Ftheta'];
BC.pdofs = BC_fixed;
nnodes = size(XYZ,1);
neq = nnodes*2;
alldofs = 1:neq;
BC.fdofs = setdiff(alldofs,BC_fixed);
BC.poi = poi;
end