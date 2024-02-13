function OPT = Optimizer( rho, FCond, BC, OPTPARAM, DOMPARAM, MATPROP,TEMPORAL_PARAM, ndof, LE, XYZ)
neq = DOMPARAM.neq;

[~, C0, ~, ~] = FEAsmallstrain(rho, FCond, BC, OPTPARAM, DOMPARAM,MATPROP,  ndof, LE, XYZ);



%% evaluate compliance
nel = DOMPARAM.nel;

% optimization PARAMETERS
x = rho;
m = 1;
n = nel;
xval=x;
epsimin = 1e-7;
xold1   = xval;
xold2   = xval;
xmin    = 0.0*ones(n, 1);
xmax    = ones(n, 1);
low     = xmin;
upp     = xmax;
c_mma   = 1000*ones(m, 1);
d_mma   = 1*ones(m, 1);
a0_mma  = 1;
a_mma   = 0*ones(m, 1);
outeriter = 0;
maxoutit  = 100;

kkttol  = 1e-3;
kktnorm = kkttol + 1;


outeriter = 0;

rmin =  OPTPARAM.rmin;
W_filt = density_filter(XYZ,LE,rmin,nel);

%% get material property parameters
E0 = MATPROP.E0;


volfrac = OPTPARAM.volfrac;


%[disp,force_val,temp,time,dfdrho]=cyclic_loading(DOMPARAM,TEMPORAL_PARAM,x,OPTPARAM,BC, 1); % IF_SENS_FLAG tells whether user wants sensitivity computed or not
while  kktnorm > kkttol %outeriter < maxoutit


    outeriter = outeriter+1;
    x_tilde=W_filt*x;

    % COMPLAINCE CONSTRAINT
    [~, C, dfdrho, KT] = FEAsmallstrain(x_tilde, FCond, BC, OPTPARAM, DOMPARAM,MATPROP, ndof, LE, XYZ);

    dCdrho = sens(x_tilde, FCond, BC, OPTPARAM, DOMPARAM, ndof, LE, XYZ, dfdrho, KT);
    dCdrho = W_filt*dCdrho;



    M=sum(x);
    M_total=DOMPARAM.nel;
    


    % OBJECTIVE FUNCTION
    f0val= C/C0;
    df0dx= dCdrho/C0;


    % VOLUME CONSTRAINT
    
    fval(1)=  M/M_total - volfrac ;
    dfdx(:,1)= (1/M_total)*ones(nel,1);


    % STIFFNESS CONSTRAINT
   
    %Eavg = sum(x.*E0)/nel;
    %fval(2)= Eavg - 0.8*E0;
    %dfdx(:,2)= (1/nel)*E0*ones(nel,1);



    % MMA FUNCTION

    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,outeriter,x,xmin,xmax,xold1,xold2,f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma);
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = x;
    x = xmma;


    %%%% The residual vector of the KKT conditions is calculated:

    [residu,kktnorm,residumax] = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);


    rhoVec = W_filt*x;
    %rhoMat = reshape(rhoVec, nely, nelx);
    PlotDensity(XYZ,LE, [], [], rhoVec, zeros(DOMPARAM.neq,1), flag)
    % figure(2)
    % imagesc(-rhoMat )
    % colormap(gray)
    % axis equal

    disp([' ************************************* ']);
    disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ...
        ' f_snap: ' sprintf('%10.7e', f0val) ' VolFrac: ' sprintf('%5.3f', M/M_total)]);
    disp([' ************************************* ']);

    density_values(outeriter,:)=rhoVec;
    save('rhoblock.mat','rhoVec');

    save('density_values.mat','density_values');
end
OPT.x = x;
end