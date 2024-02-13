


%
%
%
function main_call
clear all
close all
clc

global BC

%% XYZ matrix
XYZ_txt = readmatrix('XYZ_blisk_rev.csv');
XYZ = XYZ_txt(:,2:3);

% element conenctivity
LE_txt = readmatrix('LE_blisk_rev.csv');
LE = LE_txt(:,2:end);


% Boundary Conditions
BC_TOP = [1:12];
BC_BOTTOM = [87:101];
BC = [BC_TOP BC_BOTTOM];





% structural parameters
nnodes=size(XYZ,1);
nel=size(LE,1);


NDOF=2;
neq=nnodes*2;
penal = 1;
volfrac = 0.3;

% Plotting the airfoil design domain
mes=1;den=0;U=zeros(neq,1);
domainTRI(XYZ,LE,[],U,mes,den,[],[],[],[],penal)

% Specify the deflection of end_node required

U = zeros(neq,1);
% Calculate Area of each element

for e=1:nel
    elXY=XYZ(LE(e,:),:); % elemental nodal coordinates
    X=elXY(:,1); % elemnt X coordinate
    Y=elXY(:,2); % element Y coordinate
    a=sqrt((Y(1)-Y(2))^2+(X(1)-X(2))^2);
    b=sqrt((Y(2)-Y(3))^2+(X(2)-X(3))^2);
    c=sqrt((Y(1)-Y(3))^2+(X(1)-X(3))^2);
    s=(a+b+c)/2;
    AR(e,1)=a*b*c/(8*(s-a)*(s-b)*(s-c));
    Area(e,1)=X(1)*Y(2)+X(2)*Y(3)+X(3)*Y(1)-X(1)*Y(3)-X(2)*Y(1)-X(3)*Y(2);
end

%% 51 tail end, 37 mid top and 69 bottom
end_node=51*2;%51*2; % for airfoil_final
poi = end_node;

%%%%% define default colormap
% the boundary conditions
fixeddofs = sort([BC.*2-1; BC.*2]);
alldofs = 1:neq;
freedofs = setdiff(alldofs, fixeddofs);


% density distribution
x(1:nel,1)=volfrac;
rho = x;
% rmin calculation
rmin=3*s; %

W_filt = density_filterTRI(XYZ,LE,rmin,nel);



% Plotting the airfoil design domain
%mes=1;den=0;
%domainTRI(XYZ,LE,[],U,mes,den,[],[],[],[],penal)

DOMPARAM.neq = neq;
DOMPARAM.nel = nel;
DOMPARAM.NDOF = NDOF;
DOMPARAM.nnodes = nnodes;

BCCond.fdofs = freedofs;
BCCond.pdofs = fixeddofs;
BCCond.poi = poi;
BCCond.end_node = end_node;

OPTPARAM.penal = penal;
OPTPARAM.volfrac = volfrac;
OPTPARAM.x = x;

MATPROP.E0 = 200000; % 200, 000 MPa
MATPROP.Emin = 0;

%% evaluate the fea
Fext = zeros(neq,1);
Fext(end_node,1) = 10;

x_tilde = W_filt*x;



OPT = optimizerTRI(x_tilde, XYZ, LE, Fext, DOMPARAM,OPTPARAM,MATPROP, BCCond, Area, W_filt);
end
























%%%%%%%%%%%%%%%%%%%%%%%%%%% APPENDIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [w,x]=lgwt(N);
a=-1;
b=1;
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps


    L(:,1)=1;
    Lp(:,1)=0;

    L(:,2)=y;
    Lp(:,2)=1;

    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end

    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);

    y0=y;
    y=y0-L(:,N2)./Lp;

end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end








% function [W]=density_filter(NE,R,XYZ,LE)
% %
%
% if R == 0
%     W = eye(NE);
%     return
% end
% % %
% W = zeros(NE);
% w = zeros(NE);
% %
% for eli = 1:NE
%     ci = [(XYZ(LE(eli,1),1)+XYZ(LE(eli,3),1))/2, (XYZ(LE(eli,1),2)+XYZ(LE(eli,3),2))/2];
%     for elj = 1:NE
%         cj = [(XYZ(LE(elj,1),1)+XYZ(LE(elj,3),1))/2, (XYZ(LE(elj,1),2)+XYZ(LE(elj,3),2))/2];
%
%         r = sqrt((ci(1)-cj(1))^2 + (ci(2)-cj(2))^2);
%
%         w(eli,elj) = max(0,R-r);
%     end
%     w_norm = 1/sum(w(eli,:));
%     W(eli,:) = w_norm*w(eli,:);
%     a= sum(W(eli,:));
% end
% end


%---------------------------------------------------------------------
%  This is the file kktcheck.m
%  Version Dec 2006.
%  Krister Svanberg <krille@math.kth.se>
%
function[residu,residunorm,residumax] = ...
    kktcheck(m,n,x,y,z,lam,xsi,eta,mu,zet,s, ...
    xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
%
%  The left hand sides of the KKT conditions for the following
%  nonlinear programming problem are calculated.
%
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmax_j <= x_j <= xmin_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%   x    = Current values of the n variables x_j.
%   y    = Current values of the m variables y_i.
%   z    = Current value of the single variable z.
%  lam   = Lagrange multipliers for the m general constraints.
%  xsi   = Lagrange multipliers for the n constraints xmin_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - xmax_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general constraints.
%  xmin  = Lower bounds for the variables x_j.
%  xmax  = Upper bounds for the variables x_j.
%  df0dx = Vector with the derivatives of the objective function f_0
%          with respect to the variables x_j, calculated at x.
%  fval  = Vector with the values of the constraint functions f_i,
%          calculated at x.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at x.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%   a0   = The constants a_0 in the term a_0*z.
%   a    = Vector with the constants a_i in the terms a_i*z.
%   c    = Vector with the constants c_i in the terms c_i*y_i.
%   d    = Vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%
%*** OUTPUT:
%
% residu     = the residual vector for the KKT conditions.
% residunorm = sqrt(residu'*residu).
% residumax  = max(abs(residu)).
%
rex   = df0dx + dfdx'*lam - xsi + eta;
rey   = c + d.*y - mu - lam;
rez   = a0 - zet - a'*lam;
relam = fval - a*z - y + s;
rexsi = xsi.*(x-xmin);
reeta = eta.*(xmax-x);
remu  = mu.*y;
rezet = zet*z;
res   = lam.*s;
%
residu1 = [rex' rey' rez]';
residu2 = [relam' rexsi' reeta' remu' rezet res']';
residu = [residu1' residu2']';
residunorm = sqrt(residu'*residu);
residumax = max(abs(residu));
%---------------------------------------------------------------------
end

%-------------------------------------------------------
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
    mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
%
%    Version September 2007 (and a small change August 2008)
%
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics, SE-10044 Stockholm, Sweden.
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
%epsimin = sqrt(m+n)*10^(-9);
epsimin = 10^(-7);
raa0 = 0.00001;
albefa = 0.1;
asyinit = 0.5;
asyincr = 1.2;
asydecr = 0.7;
eeen = ones(n,1);
eeem = ones(m,1);
zeron = zeros(n,1);

% Calculation of the asymptotes low and upp :
if iter < 2.5
    low = xval - asyinit*(xmax-xmin);
    upp = xval + asyinit*(xmax-xmin);
else
    zzz = (xval-xold1).*(xold1-xold2);
    factor = eeen;
    factor(find(zzz > 0)) = asyincr;
    factor(find(zzz < 0)) = asydecr;
    low = xval - factor.*(xold1 - low);
    upp = xval + factor.*(upp - xold1);
    lowmin = xval - 10*(xmax-xmin);
    lowmax = xval - 0.01*(xmax-xmin);
    uppmin = xval + 0.01*(xmax-xmin);
    uppmax = xval + 10*(xmax-xmin);
    low = max(low,lowmin);
    low = min(low,lowmax);
    upp = min(upp,uppmax);
    upp = max(upp,uppmin);
end

% Calculation of the bounds alfa and beta :

zzz = low + albefa*(xval-low);
alfa = max(zzz,xmin);
zzz = upp - albefa*(upp-xval);
beta = min(zzz,xmax);

% Calculations of p0, q0, P, Q and b.

xmami = xmax-xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
xmamiinv = eeen./xmami;
ux1 = upp-xval;
ux2 = ux1.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
%
p0 = zeron;
q0 = zeron;
p0 = max(df0dx,0);
q0 = max(-df0dx,0);
%p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
%q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
p0 = p0 + pq0;
q0 = q0 + pq0;
p0 = p0.*ux2;
q0 = q0.*xl2;
%
P = sparse(m,n);
Q = sparse(m,n);
P = max(dfdx,0);
Q = max(-dfdx,0);
%P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
%Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
P = P + PQ;
Q = Q + PQ;
P = P * spdiags(ux2,0,n,n);
Q = Q * spdiags(xl2,0,n,n);
b = P*uxinv + Q*xlinv - fval ;
%
%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
    subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);




end

%-------------------------------------------------------------
%    This is the file subsolv.m
%
%    Version Dec 2006.
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics, KTH,
%    SE-10044 Stockholm, Sweden.
%
function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
    subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
%
% This function subsolv solves the MMA subproblem:
%
% minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
%          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
%
% subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
%            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
%
% Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
% Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
%
een = ones(n,1);
eem = ones(m,1);
epsi = 1;
epsvecn = epsi*een;
epsvecm = epsi*eem;
x = 0.5*(alfa+beta);
y = eem;
z = 1;
lam = eem;
xsi = een./(x-alfa);
xsi = max(xsi,een);
eta = een./(beta-x);
eta = max(eta,een);
mu  = max(eem,0.5*c);
zet = 1;
s = eem;
itera = 0;
while epsi > epsimin
    epsvecn = epsi*een;
    epsvecm = epsi*eem;
    ux1 = upp-x;
    xl1 = x-low;
    ux2 = ux1.*ux1;
    xl2 = xl1.*xl1;
    uxinv1 = een./ux1;
    xlinv1 = een./xl1;
    plam = p0 + P'*lam ;
    qlam = q0 + Q'*lam ;
    gvec = P*uxinv1 + Q*xlinv1;
    dpsidx = plam./ux2 - qlam./xl2 ;
    rex = dpsidx - xsi + eta;
    rey = c + d.*y - mu - lam;
    rez = a0 - zet - a'*lam;
    relam = gvec - a*z - y + s - b;
    rexsi = xsi.*(x-alfa) - epsvecn;
    reeta = eta.*(beta-x) - epsvecn;
    remu = mu.*y - epsvecm;
    rezet = zet*z - epsi;
    res = lam.*s - epsvecm;
    residu1 = [rex' rey' rez]';
    residu2 = [relam' rexsi' reeta' remu' rezet res']';
    residu = [residu1' residu2']';
    residunorm = sqrt(residu'*residu);
    residumax = max(abs(residu));
    ittt = 0;
    while residumax > 0.9*epsi & ittt < 200
        ittt=ittt + 1;
        itera=itera + 1;
        ux1 = upp-x;
        xl1 = x-low;
        ux2 = ux1.*ux1;
        xl2 = xl1.*xl1;
        ux3 = ux1.*ux2;
        xl3 = xl1.*xl2;
        uxinv1 = een./ux1;
        xlinv1 = een./xl1;
        uxinv2 = een./ux2;
        xlinv2 = een./xl2;
        plam = p0 + P'*lam ;
        qlam = q0 + Q'*lam ;
        gvec = P*uxinv1 + Q*xlinv1;
        GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
        dpsidx = plam./ux2 - qlam./xl2 ;
        delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
        dely = c + d.*y - lam - epsvecm./y;
        delz = a0 - a'*lam - epsi/z;
        dellam = gvec - a*z - y - b + epsvecm./lam;
        diagx = plam./ux3 + qlam./xl3;
        diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
        diagxinv = een./diagx;
        diagy = d + mu./y;
        diagyinv = eem./diagy;
        diaglam = s./lam;
        diaglamyi = diaglam+diagyinv;
        if m < n
            blam = dellam + dely./diagy - GG*(delx./diagx);
            bb = [blam' delz]';
            Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
            AA = [Alam     a
                a'    -zet/z ];
            solut = AA\bb;
            dlam = solut(1:m);
            dz = solut(m+1);
            dx = -delx./diagx - (GG'*dlam)./diagx;
        else
            diaglamyiinv = eem./diaglamyi;
            dellamyi = dellam + dely./diagy;
            Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
            azz = zet/z + a'*(a./diaglamyi);
            axz = -GG'*(a./diaglamyi);
            bx = delx + GG'*(dellamyi./diaglamyi);
            bz  = delz - a'*(dellamyi./diaglamyi);
            AA = [Axx   axz
                axz'  azz ];
            bb = [-bx' -bz]';
            solut = AA\bb;
            dx  = solut(1:n);
            dz = solut(n+1);
            dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        end
        %
        dy = -dely./diagy + dlam./diagy;
        dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
        deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
        dmu  = -mu + epsvecm./y - (mu.*dy)./y;
        dzet = -zet + epsi/z - zet*dz/z;
        ds   = -s + epsvecm./lam - (s.*dlam)./lam;
        xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
        dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
        %
        stepxx = -1.01*dxx./xx;
        stmxx  = max(stepxx);
        stepalfa = -1.01*dx./(x-alfa);
        stmalfa = max(stepalfa);
        stepbeta = 1.01*dx./(beta-x);
        stmbeta = max(stepbeta);
        stmalbe  = max(stmalfa,stmbeta);
        stmalbexx = max(stmalbe,stmxx);
        stminv = max(stmalbexx,1);
        steg = 1/stminv;
        %
        xold   =   x;
        yold   =   y;
        zold   =   z;
        lamold =  lam;
        xsiold =  xsi;
        etaold =  eta;
        muold  =  mu;
        zetold =  zet;
        sold   =   s;
        %
        itto = 0;
        resinew = 2*residunorm;
        while resinew > residunorm & itto < 50
            itto = itto+1;
            x   =   xold + steg*dx;
            y   =   yold + steg*dy;
            z   =   zold + steg*dz;
            lam = lamold + steg*dlam;
            xsi = xsiold + steg*dxsi;
            eta = etaold + steg*deta;
            mu  = muold  + steg*dmu;
            zet = zetold + steg*dzet;
            s   =   sold + steg*ds;
            ux1 = upp-x;
            xl1 = x-low;
            ux2 = ux1.*ux1;
            xl2 = xl1.*xl1;
            uxinv1 = een./ux1;
            xlinv1 = een./xl1;
            plam = p0 + P'*lam ;
            qlam = q0 + Q'*lam ;
            gvec = P*uxinv1 + Q*xlinv1;
            dpsidx = plam./ux2 - qlam./xl2 ;
            rex = dpsidx - xsi + eta;
            rey = c + d.*y - mu - lam;
            rez = a0 - zet - a'*lam;
            relam = gvec - a*z - y + s - b;
            rexsi = xsi.*(x-alfa) - epsvecn;
            reeta = eta.*(beta-x) - epsvecn;
            remu = mu.*y - epsvecm;
            rezet = zet*z - epsi;
            res = lam.*s - epsvecm;
            residu1 = [rex' rey' rez]';
            residu2 = [relam' rexsi' reeta' remu' rezet res']';
            residu = [residu1' residu2']';
            resinew = sqrt(residu'*residu);
            steg = steg/2;
        end
        residunorm=resinew;
        residumax = max(abs(residu));
        steg = 2*steg;
    end
    if ittt > 198
        epsi
        ittt
    end
    epsi = 0.1*epsi;
end
xmma   =   x;
ymma   =   y;
zmma   =   z;
lamma =  lam;
xsimma =  xsi;
etamma =  eta;
mumma  =  mu;
zetmma =  zet;
smma   =   s;
%-------------------------------------------------------------

end










