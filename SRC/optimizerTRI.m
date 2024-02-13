function OPT = optimizerTRI(x_tilde, XYZ, LE, Fext, DOMPARAM,OPTPARAM,MATPROP, BCCond, Area, W_filt)
nel = DOMPARAM.nel;
neq = DOMPARAM.neq;

[Fint_d1, Ktan, Fint, U] = FEATri(x_tilde, XYZ, LE, Fext, DOMPARAM,OPTPARAM,MATPROP, BCCond);

%% evaluate the compliance and the sensitivities
[C0, dCdrho]=complianceTRI(XYZ, LE, U, Fext, Ktan, Fint_d1, DOMPARAM, BCCond);


% calculate constraint function
outeriter=0;
kkttol  = 1e-5;
kktnorm = kkttol + 1;
kktcond=kktnorm;
volfrac = OPTPARAM.volfrac;
x = OPTPARAM.x;
end_node = BCCond.end_node;
% MASS CONSTRAINT
M_bound=volfrac;
M=sum(diag(Area)*(x));
M_total=sum(Area);
vol=M/M_total;



% optimization PARAMETERS
xval=x;
m = 1;
n=nel;
epsimin = 1e-7;
xold1   = xval;
xold2   = xval;
xmin    = 0.1*ones(n, 1);
xmax    = ones(n, 1);
low     = xmin;
upp     = xmax;
c_mma   = 1000*ones(m, 1);
d_mma   = 1*ones(m, 1);
a0_mma  = 1;
a_mma   = 0*ones(m, 1);

maxoutit  = 100;

Elim = 0.8*(1/nel)*sum(ones(nel,1).*MATPROP.E0);
while  kktnorm > kkttol %outeriter < maxoutit



    
    rhoVec =W_filt*x;
    U_end=U(end_node);
    

    outeriter = outeriter+1;


    x_tilde=W_filt*x;



    x_tilde = W_filt*x;
    [Fint_d1, Ktan, Fint, U] = FEATri(x_tilde, XYZ, LE, Fext, DOMPARAM,OPTPARAM,MATPROP, BCCond);

    %% evaluate the compliance and the sensitivities
    [C, dCdrho]=complianceTRI(XYZ, LE, U, Fext, Ktan, Fint_d1, DOMPARAM, BCCond);
    dCdrho = W_filt*dCdrho;



    M=sum(diag(Area)*(x));
    M_total=sum(Area);
    vol=M/M_total;
    %   vol=vol/M_bound

    % OBJECTIVE FUNCTION
    f0val= C/C0;
    df0dx= dCdrho/C0;

    

    % VOLUME CONSTRAINT
%     fval(1)= M/M_total - M_bound;
%     dfdx(:,1)= (1/M_total)*Area;
  
    %% Elasticity Constraint
    Eavg = sum(x.*MATPROP.E0)/nel;
    fval(1)= Eavg - Elim;
    dfdx(:,1)= (1/nel)*(OPTPARAM.penal*x.^(OPTPARAM.penal-1))*(MATPROP.E0);
    
    

    % MMA FUNCTION

    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,outeriter,x,xmin,xmax,xold1,xold2,f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma);
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = x;
    x = xmma;


    %%%% The residual vector of the KKT conditions is calculated:

    [residu,kktnorm,residumax] = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);

% Plot material distribution
    disp([' ************************************* ']);
    disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ...
        ' f_snap: ' sprintf('%10.7e', f0val) ' VolFrac: ' sprintf('%5.3f', M/M_total) ' ' num2str(U(end_node))]);
    disp([' ************************************* ']);


    % %     x_tilde=W_filt*x;
    mes=0;den=1;U=zeros(neq,1);
    domainTRI(XYZ,LE,rhoVec,U,mes,den,outeriter,kktnorm,vol,U_end,OPTPARAM.penal)


    density_values(outeriter,:)=rhoVec;
    save('rhoblock.mat','rhoVec');

    save('density_values.mat','density_values');
end
OPT.x = x;
end