%%%%%%%%%%%%%%%   TOTAL STIFFNESS MATRIX      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Fint_d1,Ktan,Fint]= stiffTRI(nel,nnodes,NDOF,LE,XYZ,u,x,penal, MATPROP)
neq=nnodes*NDOF;
Kgeo   = speye(neq,neq);
Kmat   = speye(neq,neq);
J      = sparse(neq);
Fint   = zeros(neq,1);
Fint_d1=zeros(neq,nel);

for e=1:nel

    [fint_d1,fe, kmat] = ElementEvalTri(e,u,XYZ,LE,nel,x,penal, MATPROP);
    % output=struct('array1',fint_d1,'array2',fe,'array3',kgeo,'array4',kmat);


    for I=1:3
        II=(I-1)*NDOF+1;
        IDOF(II:II+1)=(LE(e,I)-1)*NDOF+1:(LE(e,I)-1)*NDOF+2;
    end


    Fint(IDOF)      = Fint(IDOF)+fe;
    Fint_d1(IDOF,e)   = Fint_d1(IDOF,e)+ fint_d1;
    %
    
    Kmat(IDOF,IDOF) = Kmat(IDOF,IDOF) + kmat;
    %      J(IDOF,IDOF)    = J(IDOF,IDOF) + je;
end


Kmat = Kmat- speye(neq);


Ktan = Kmat ;
end