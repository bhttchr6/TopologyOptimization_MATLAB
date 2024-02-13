function dCdrho = sens(rho, FCond, BC, OPTPARAM, DOMPARAM, ndof, LE, XYZ, Fint_d1, KT)
poi = BC.poi;
neq = DOMPARAM.neq;
nel = DOMPARAM.nel;
lambda = zeros(neq,1);
L = zeros(neq,1);
L(poi,1) = 1;
L_fT = L(BC.fdofs,1)';
%Fext = ForceEval(FCond, DOMPARAM.neq);
Fext = ApplyForce(FCond, DOMPARAM, "along y");
Fext_f = Fext(BC.fdofs,1)';





kff = KT(BC.fdofs, BC.fdofs);
kpf = KT(BC.pdofs, BC.fdofs);
lambdaf_T= -Fext_f/(kff);
%lambdaf_T= -L_fT/(kff);

%
lambda(BC.fdofs,1)= lambdaf_T';


for IE=1:nel

    IDOF=zeros(1,8);
    for I=1:4
        II=(I-1)*ndof+1;
        IDOF(II:II+1)=(LE(IE,I)-1)*ndof+1:(LE(IE,I)-1)*ndof+2;
    end

    lambda_e=lambda(IDOF);

    dc(IE,:)=lambda_e'*Fint_d1(IDOF,IE);
end
dCdrho = dc;
% %
% rho_orig = rho;
% h = 1e-06;
% [~, C_0, ~, ~] = FEAsmallstrain(rho_orig, FCond, BC, OPTPARAM, DOMPARAM, ndof, LE, XYZ);
% 
% for e = 1:size(rho,1)
%     rho(e,1) = rho(e,1)+ h;
%     [~, C_p, ~, ~] = FEAsmallstrain(rho, FCond, BC, OPTPARAM, DOMPARAM, ndof, LE, XYZ);
%     dCdrho_FD(e,1) = (C_p - C_0)/h;
%     rho(e,1) = rho_orig(e,1);
% end
% 
% diff = max(dc-dCdrho_FD);

end