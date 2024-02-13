function [Fint_d1,Ktan,Fint]= stiff(rho, OPTPARAM, DOMPARAM,MATPROP, NDOF,LE,XYZ,u)

neq = DOMPARAM.neq;
nel = DOMPARAM.nel;

Kmat   = speye(neq,neq);
Fint   = zeros(neq,1);
Fint_d1=zeros(neq,nel);


for e=1:nel
    
    [fint_d1,fe, kmat] = ElementEval(rho,OPTPARAM,MATPROP, e,u,XYZ,LE,nel);
    
    for I=1:4
        II=(I-1)*NDOF+1;
        IDOF(II:II+1)=(LE(e,I)-1)*NDOF+1:(LE(e,I)-1)*NDOF+2;
    end
    Fint(IDOF)      = Fint(IDOF) + fe;
    Fint_d1(IDOF,e)   = Fint_d1(IDOF,e)+ fint_d1;
    
    
    
    Kmat(IDOF,IDOF) = Kmat(IDOF,IDOF) + kmat;
   
end

Kmat = Kmat- speye(neq);

Ktan = Kmat;
end