%%% Small deformation TO code
function main
clear all
close all
clc

%% Define the degrees of freedom
ndof = 2;


%% Read GMSH based mesh (save mesh as .inp and then separate nodes/elements)
XYZ = readmatrix("GEOM/TXT files/XYZ_impact.csv");
XYZ = XYZ(:, 2:end-1);

LE = readmatrix("GEOM/TXT files/LE_impact.csv");
LE = LE(:, 2:end);

%%% define XYZ
%[XYZ, LE, DOMPARAM] = Domain(nelx, nely, L, H, ndof);
neq = size(XYZ,1)*ndof;
DOMPARAM.neq = neq;
DOMPARAM.nel = size(LE,1);
DOMPARAM.eps = 0.5;
%PlotBVP(XYZ,LE, [], [], [], zeros(DOMPARAM.neq,1), 1);

%% Define BC
bottom_BC = readmatrix("GEOM/TXT files/BC.csv");
bottom_BC = sort(bottom_BC(~isnan(bottom_BC)));

right_BC = readmatrix("GEOM/TXT files/right_BC.csv");
right_BC = sort(right_BC(~isnan(right_BC)));

right_peak_BC = readmatrix("GEOM/TXT files/right_peak_BC.xlsx");
right_peak_BC = sort(right_peak_BC(~isnan(right_peak_BC)));

valley_BC = readmatrix("GEOM/TXT files/valley_BC.csv");
valley_BC = sort(valley_BC(~isnan(valley_BC)));

left_peak_BC = readmatrix("GEOM/TXT files/left_peak_BC.xlsx");
left_peak_BC = sort(left_peak_BC(~isnan(left_peak_BC)));

left_BC = readmatrix("GEOM/TXT files/left_BC.csv");
left_BC = sort(left_BC(~isnan(left_BC)));

pdofs1 = ApplyBC(right_peak_BC, "fixed y");
pdofs2 = ApplyBC(right_peak_BC, "fixed x");

pdofs3 = ApplyBC(left_peak_BC, "fixed y");
pdofs4 = ApplyBC(left_peak_BC, "fixed x");

pdofs5 = ApplyBC(left_BC, "roller y");
pdofs6 = ApplyBC(right_BC, "roller y");

BC_fixed_dofs = [pdofs1' pdofs2' pdofs3' pdofs4' pdofs5' pdofs6'];
alldofs = 1:neq;
freedofs = setdiff(alldofs, BC_fixed_dofs);

BC.pdofs = BC_fixed_dofs;
BC.fdofs = freedofs;
BC.poi = valley_BC.*2;


%% Force Parameters
ForceMagnitude = 1;
FCond.ForceMagnitude = ForceMagnitude;
FCond.ForceNodes = bottom_BC;
Force =  ApplyForce(FCond, DOMPARAM, "along y");
%% define BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BVP =1, 101 (cantilever)
% BVP = 2(MBB beam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BVP = 2;

%theta_0 = 270;
%[FCond, BC_fixed_dofs, BC] = BoundaryConditions(nelx, nely, XYZ, LE, ndof, BVP, ForceMagnitude, theta_0);

%% visualize domain
%PlotBVP(XYZ,LE, BC_fixed_dofs,bottom_BC , [], zeros(DOMPARAM.neq,1), 1);


%% Optimization
rho = ones(DOMPARAM.nel,1);
OPTPARAM.penal = 3;
OPTPARAM.Emin = 1e-06;
OPTPARAM.rmin = 2*DOMPARAM.eps;
OPTPARAM.volfrac = 0.3;

MATPROP.E0 = 200000;
%PlotBVP(XYZ,LE, BC_fixed_dofs, FCond(:,1), FCond(:,2), U, 2);

%% Set the temporal specifications
TEMPORAL_PARAM.time_spec = 35;
TEMPORAL_PARAM.delta_t = 5;

OPT = Optimizer(rho, FCond, BC, OPTPARAM, DOMPARAM, MATPROP, TEMPORAL_PARAM, ndof, LE, XYZ);

%[U, C0, ~, ~] = FEAsmallstrain(rho, FCond, BC, OPTPARAM, DOMPARAM,MATPROP,  ndof, LE, XYZ);

%PlotBVP(XYZ,LE, BC_fixed_dofs,bottom_BC , [], U, 1);

end
