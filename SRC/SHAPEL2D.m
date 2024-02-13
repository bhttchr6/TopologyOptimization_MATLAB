%function [N,B, DET,dN_dx] = SHAPEL2D(psi,eta,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x=C(:,1);
% y=C(:,2);
% N=[0.25*(1-psi)*(1-eta) 0 0.25*(1+psi)*(1-eta) 0  0.25*(1+psi)*(1+eta) 0 0.25*(1-psi)*(1+eta) 0;
%     0 0.25*(1-psi)*(1-eta) 0 0.25*(1+psi)*(1-eta) 0 0.25*(1+psi)*(1+eta) 0 0.25*(1-psi)*(1+eta)];
% % %
% dN_dshi=[-0.25*(1-eta) 0.25*(1-eta) 0.25*(1+eta) -0.25*(1+eta)];
% dN_deta=[-0.25*(1-psi) -0.25*(1+psi) 0.25*(1+psi) 0.25*(1-psi)];
% 
% %
% dx_dshi=dN_dshi*x;
% dx_deta=dN_deta*x;
% dy_dshi=dN_dshi*y;
% dy_deta=dN_deta*y;
% Jacobian=[dx_dshi dy_dshi;dx_deta dy_deta];
% 
% dN_dx=inv(Jacobian)*[dN_dshi;dN_deta];
% DET=det(Jacobian);
% 
% B=[dN_dx(1,1) 0 dN_dx(1,2) 0 dN_dx(1,3) 0 dN_dx(1,4) 0;
%     0       dN_dx(2,1) 0 dN_dx(2,2) 0 dN_dx(2,3) 0 dN_dx(2,4);
%     dN_dx(2,1) dN_dx(1,1) dN_dx(2,2) dN_dx(1,2) dN_dx(2,3) dN_dx(1,3) dN_dx(2,4) dN_dx(1,4)];

function [B, DET] = SHAPEL2D(psi,eta,C)

GN    = 0.25 * [eta-1  1-eta   1+eta   -eta-1;
    psi-1  -psi-1  1+psi    1-psi];
J     = GN*C;        % Get the Jacobian matrix

DET  = det(J);     % Needed for the coordinate transformations

BB     = J\GN;       % compute the derivative of the shape functions

B1x     = BB(1,1);
B2x     = BB(1,2);
B3x     = BB(1,3);
B4x     = BB(1,4);
B1y     = BB(2,1);
B2y     = BB(2,2);
B3y     = BB(2,3);
B4y     = BB(2,4);



B = [ B1x      0     B2x     0      B3x    0      B4x     0  ;
    0     B1y     0     B2y      0     B3y     0      B4y;
    B1y     B1x    B2y    B2x     B3y    B3x    B4y     B4x];





end