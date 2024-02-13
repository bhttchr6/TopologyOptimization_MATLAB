function [fint_d1,fint,kmat]= ElementEvalTri(e,u,XYZ,LE,nel,x,penal, MATPROP)

pr = 0.3;
ngp=1; % no. of gauss points in each direction
nen=3;
NDOF=2;
IE=e;

E0 = MATPROP.E0;
Emin = MATPROP.Emin;%1e-06;
E = x(e,1)^penal*E0 + Emin;
ndof=NDOF;
dEdrho = penal*x(e,1)^(penal-1)*E0;

ind2voit=[1; 4; 3];

elXY=XYZ(LE(e,:),:); % elemental nodal coordinates
X=elXY(:,1); % elemnt X coordinate
Y=elXY(:,2); % element Y coordinate

IDOF=zeros(1,6);
for I=1:3
    II=(I-1)*NDOF+1;
    IDOF(II:II+1)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+2;
end
ue=u(IDOF); % elemnt nodal displacements
ux=[ue(1) ue(3) ue(5)]; % X displacements
uy=[ue(2) ue(4) ue(6)]; % Y displacements

Ccur=[X+ux' Y+uy']; % current coordinates

[w,gp] = lgwt(ngp); % weights and integration points

kgeo=zeros(nen*ndof,nen*ndof); % stiffness matrix
kmat=zeros(nen*ndof,nen*ndof);
fint=zeros(nen*ndof,1);        % internal force
fint_d1=zeros(nen*ndof,1);


for i=1:ngp
    for j=1:ngp
        eta = gp(i);
        psi = gp(j);

        [dNdx,Bx, detJx]     = SHAPEL2DTRI(Ccur);     % derivative of the shape functions

        strain = Bx*ue;

        
        %% Stiffness matrix
        D = E/((1+pr)*(1-2*pr))*[1-pr pr   0;
                                  pr  1-pr 0;
                                  0   0    0.5-pr];

        
        dDdrho = dEdrho/((1+pr)*(1-2*pr))*[1-pr pr   0;
                                            pr  1-pr 0;
                                            0   0    0.5-pr];

        cauchy = D*strain;
        dcauchydrho = dDdrho*strain;
        fint = fint + Bx'*cauchy*detJx*w(i)*w(j);
        fint_d1=fint_d1+Bx'*dcauchydrho*detJx*w(i)*w(j);

 
        kmat = kmat + w(i)*w(j)*Bx'*D*Bx*detJx;   % material element matrix



    end
end

end