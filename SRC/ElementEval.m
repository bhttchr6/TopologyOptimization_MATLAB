function [fint_d1,fint,kmat]=ElementEval(rho, OPTPARAM,MATPROP, e,u,XYZ,LE,nel)

ind2voit=[1; 4; 3];
penal = OPTPARAM.penal;
Emin = OPTPARAM.Emin;


pr = 0.3;
ngp = 2; % no. of gauss points in each direction
nen = 4;
NDOF = 2;
IE=e;

ndof=NDOF;

E0 = MATPROP.E0; %% 200, 000 Mpa
E = rho(e,1)^penal*E0 + Emin;
dEdrho = penal*rho(e,1)^(penal-1)*E0;

elXY=XYZ(LE(e,:),:); % elemental nodal coordinates
X=elXY(:,1); % elemnt X coordinate
Y=elXY(:,2); % element Y coordinate

IDOF=zeros(1,8);
for I=1:4
    II=(I-1)*NDOF+1;
    IDOF(II:II+1)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+2;
end

ue=u(IDOF); % elemnt nodal displacements
ux=[ue(1) ue(3) ue(5) ue(7)]; % X displacements
uy=[ue(2) ue(4) ue(6) ue(8)]; % Y displacements

Ccur=[X Y]; % current coordinates



[w,gp] = lgwt(ngp); % weights and integration points


kmat=zeros(nen*ndof,nen*ndof);
fint=zeros(nen*ndof,1);        % internal force
fint_d1=zeros(nen*ndof,1);

for i=1:ngp
    for j=1:ngp
        eta = gp(i);
        psi = gp(j);

        [B, DET]     = SHAPEL2D(psi,eta,Ccur);

        strain = B*ue;
        %% Stiffness matrix
        D = E/((1+pr)*(1-2*pr))*[1-pr pr   0;
                                  pr  1-pr 0;
                                  0   0    0.5-pr];

        
        dDdrho = dEdrho/((1+pr)*(1-2*pr))*[1-pr pr   0;
                                            pr  1-pr 0;
                                            0   0    0.5-pr];

        cauchy = D*strain;
        dcauchydrho = dDdrho*strain;
        

        %cauchy_bar=[cauchy(1,1);cauchy(2,2);cauchy(1,2)];


        fint =fint+B'*cauchy*DET*w(i)*w(j);

        fint_d1 = fint_d1 + B'*dcauchydrho*DET*w(i)*w(j);

        kmat = kmat + w(i)*w(j)*B'*D*B*DET; % material element matrix



    end
end



end