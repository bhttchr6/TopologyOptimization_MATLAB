function [C, dCdrho]=complianceTRI(XYZ, LE, U, Fext, KT, Fint_d1, DOMPARAM, BCCond)

C = sum(Fext.*U);

poi = BCCond.poi;
neq = DOMPARAM.neq;
nel = DOMPARAM.nel;
NDOF = DOMPARAM.NDOF;

lambda = zeros(neq,1);
L = zeros(neq,1);
L(poi,1) = 1;
L_fT = L(BCCond.fdofs,1)';
%Fext = ForceEval(FCond, DOMPARAM.neq);
Fext_f = Fext(BCCond.fdofs,1)';





kff = KT(BCCond.fdofs, BCCond.fdofs);
kpf = KT(BCCond.pdofs, BCCond.fdofs);
lambdaf_T= -Fext_f/(kff);
%lambdaf_T= -L_fT/(kff);

%
lambda(BCCond.fdofs,1)= lambdaf_T';


for IE=1:nel

    IDOF=zeros(1,6);
    for I=1:3
        II=(I-1)*NDOF+1;
        IDOF(II:II+1)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+2;
    end

    lambda_e=lambda(IDOF);

    dc(IE,:)=lambda_e'*Fint_d1(IDOF,IE);
end
dCdrho = dc;
end