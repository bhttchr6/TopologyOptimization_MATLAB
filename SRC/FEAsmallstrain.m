function [U, C, Fint_d1, Ktan] = FEAsmallstrain(rho, FCond, BC, OPTPARAM, DOMPARAM, MATPROP, ndof, LE, XYZ)
%load_node = BC.poi;
Fext = ApplyForce(FCond, DOMPARAM, "along y");

U = zeros(DOMPARAM.neq,1);
[Fint_d1,Ktan,Fint]= stiff(rho, OPTPARAM, DOMPARAM,MATPROP, ndof,LE,XYZ,U);

residual = Fext-Fint;
U(BC.fdofs,1) = Ktan(BC.fdofs, BC.fdofs)\residual(BC.fdofs,1);

%PlotBVP(XYZ,LE, [], [], [], U, 2);
%% calculate internal force
[Fint_d1,Ktan,Fint]= stiff(rho, OPTPARAM, DOMPARAM,MATPROP, ndof,LE,XYZ,U);
%% complaince of the system
%C = sum(Fext.*U);
C = sum(U(BC.poi,1));
end