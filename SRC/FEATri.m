function [Fint_d1, Ktan, Fint, U] = FEATri(rho, XYZ, LE, Fext, DOMPARAM,OPTPARAM,MATPROP,  BCCond)

%rho = ones(nel,1);
x = rho;
U = zeros(DOMPARAM.neq, 1);
nel = DOMPARAM.nel;
nnodes = DOMPARAM.nnodes;
NDOF = DOMPARAM.NDOF;
penal = OPTPARAM.penal;

freedofs = BCCond.fdofs;
fixeddofs = BCCond.pdofs;

[ ~,Ktan,Fint]= stiffTRI(nel,nnodes,NDOF,LE,XYZ,U,x,penal, MATPROP);



Residual = Fext-Fint;



U(freedofs,1) = Ktan(freedofs,freedofs)\Residual(freedofs,1);

[ Fint_d1,Ktan,Fint]= stiffTRI(nel,nnodes,NDOF,LE,XYZ,U,x,penal, MATPROP);
end