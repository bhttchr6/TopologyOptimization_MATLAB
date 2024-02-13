% MAIN FORWARD ANALYSIS AND SESNITIVITY ANALYSIS FUNCTION
%function [disp,force_val,temp,time,dfdrho]=cyclic_loading(nelx,nely,time_spec,delta_t,x,penal,sens)
function [disp,force_val,temp,time,dfdrho]=cyclic_loading(DOMPARAM,TEMPORAL_PARAM,x,OPTPARAM,BC, IF_SENS_FLAG)

global Tmax Tmin Tg XYZ LE ndim   load_dof
%Tmax=355;
% Tmin=275;
%Tmin=325;
Tmax=350;
%Tg=343;
Tg=340;
Tmin=330;

ndim=2;








%  load_dof=5;
load_dof = BC.poi;  % for shear
% load_dof=nnodes*2-((n_nodes_y)*2-1);
% u_marker=0.2*l;  %20 percent strain
u_marker=0.2*DOMPARAM.H;  %20 percent strain % for shear





%struc_param=[nel nnodes neq n_nodes_x n_nodes_y ndof];







% initialize the SMP variables
%%%%%%%% FORWARD ANALYSIS %%%%%%%%%%%%%%

iter=0;

iter_outer_loop=0;

U=zeros(DOMPARAM.neq,1);
U_hist=zeros(DOMPARAM.neq,1);
T_hist=Tmax;
% theta=0;
% 
% % specify Nodal parameters
alldofs = 1:DOMPARAM.neq;
fixeddofs = BC.pdofs;% half MBB with fixed ends
freedofs = BC.fdofs;
nbc=find(fixeddofs<load_dof);
nbc=size(nbc,2);
% 
% 
%% Stimulus function %%%%%%%%%%%%
STIMULUS_PARAM.Tmax = Tmax;
STIMULUS_PARAM.Tmin = Tmin;
tinitial = 0;
delta_t = TEMPORAL_PARAM.delta_t;
HR = 1;CR = 1; 
STIMULUS_PARAM.tinitial = tinitial;
STIMULUS_PARAM.HR = HR;
STIMULUS_PARAM.CR = CR;
STIMULUS_PARAM.t_final = TEMPORAL_PARAM.time_spec;
n_timesteps = 50;%(STIMULUS_PARAM.t_final - STIMULUS_PARAM.tinitial)/delta_t;
t = tinitial;
for titer = 1:n_timesteps
    t = t + delta_t;
T(titer) = stimulus_function(STIMULUS_PARAM, t);
end
figure()
plot(1:n_timesteps, T);

% %%%%%%%%%%%%%%%%%%%%%%%%%%  STEP-III: COOLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 delta_t=5;
 val=0.005; 
 Fext=zeros(neq,1);
% % FORCE CONTROL METHOD
%  load('workspace','U','Fext','freedofs_I')
%  delta_t=5;
% U_hist=U;
% T_hist=Tmax;
  T=Tmax;
 C_R=-1;
 iter_cool=iter+1;
t=0;
t_cool_start=t;
%  t_cool_start=t_relaxI_fin;
t_cool_fin=t_cool_start+(Tmax-Tmin)/-C_R;
%t=0;
% t_cool_start=t;
% t_cool_fin=t_cool_start+(Tmax-Tmin)/-C_R;
if time_spec>t_cool_start
    
    while t+delta_t<=min(t_cool_fin,time_spec)
         t=t+delta_t;
        
        %         if t<=50
                      Fext(load_dof)=val;
        %         end
        %         if t>50
        %                Fext(load_dof)=0.00;
        %             Fext=zeros(neq,1);
        %         end
        %
        
        %          T=Tmax;
        T=T+C_R*delta_t;
        h_c=1;
        iter=iter+1;
        
        
        
        
        [~,~,Ktan,Khistry, Fint,~]=assembly_SMP(T,U,U_hist,delta_t,T_hist,struc_param,h_c,x,penal,iter,0);
%         Fext=Fint;
                residual=Fint-Fext;
                res=norm(residual(freedofs,:));
                % res=1;
                while res>1e-12
                    %              [~,~,Ktan,Fint]=assembly_SMP(T,U,U_hist,delta_t,T_hist,struc_param,h_c,x,penal,iter,0);
                    %              residual=Fint-Fext;
                    %             res=norm(residual(freedofs,:));
                    dU=-Ktan(freedofs,freedofs)\(residual(freedofs,:));
                    U(freedofs,:)=U(freedofs,:)+dU;
                    [~,~,Ktan,Khistry,Fint,D]=assembly_SMP(T,U,U_hist,delta_t,T_hist,struc_param,h_c,x,penal,iter,0);
                    residual=Fint-Fext;
                    res=norm(residual(freedofs,:));
                end
% %                 Fext=Fint;
%         
%                                   if iter==1
%                                       h=1e-06;
%                                       x(1,1)=x(1,1)+h;
%          
%                                   end
        %                     if iter==6
        % %                         U_hist(:,iter)
        %                          h=rand(1,1).*1e-06
        %                          U_hist(8,iter)=U_hist(8,iter)+h;
        %                            U_hist(:,iter)
        %                    end
        % calculate sensitivity of each step after residual satisfaction
        [dFint_dun,dFint_dx,~,~,Fint,~]=assembly_SMP(T,U,U_hist,delta_t,T_hist,struc_param,h_c,x,penal,iter,1);
        
        %        if iter==6
        %
        %            Fint
        %           dFint_dx
        %           dFint_dun
        %        end
        Fext=Fint;
        U_hist=[U_hist U];
        T_hist=[T_hist T];
        % update variables
        
        
        %         stress(iter)=Fint(load_dof,1);
%         t=t+delta_t;
        temp(iter)=T;
        force_val(iter)=Fext(load_dof);
        time(iter)=t;
        disp(iter)=U(load_dof);
        
        %%% store tangent stiffness values
        KTAN(:,:,iter)=Ktan;
        KHIST(:,:,iter)=Khistry;
        dFintdun(:,:,iter)=dFint_dun;
        dFintdx(:,:,iter)=dFint_dx;
    end
%                   figure(1)
%              plot(time,disp,'linewidth',2)
%              xlabel('time')
%              ylabel('displacement')
%              grid on
    
             % domain(XYZ,LE,[],U,1,[],[],[],[],[],[])
    
end




fprintf("FEA COMPLETED");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SENSITIVITY ANALYSIS
dfdrho=zeros(nel,1);
if sens==1
    iter_final=iter;
    dfdrho=zeros(nel,1);
    lambda=zeros(neq,iter);
    %
    % %%% compute final adjoint vector by firstprinciple
    %
    L=zeros(neq,1);
    L(load_dof)=1;
    p=size(fixeddofs,2);
    f=size(freedofs,2);
    dFextmdYm=zeros(neq,neq);
    dFintmdYm=zeros(neq,neq);
    dFintdY=zeros(neq,neq);
    dFextdY=zeros(neq,neq);
    dFintdYprev=zeros(neq,neq);
    
    
    for i=1:size(fixeddofs,2)
        val1=fixeddofs(i);
        dFextmdYm(val1,val1)=1;
    end
%     dFextmdYm(1:p,1:p)=eye(p);
    
    for i=1:size(freedofs,2)
        val2=freedofs(i);
        
        dFintmdYm(:,val2)=KTAN(:,val2,iter);
    end

    dRmdYm=dFintmdYm-dFextmdYm;
    lambda_M=-L'/[dRmdYm];
    dRmdYm_inv=inv(dRmdYm);
    N=nnz(dRmdYm(freedofs,freedofs));
    lambda(:,iter)=lambda_M;
%-----------------------------------------------------------------------    
%     %%% manual adjoint formulations
%     %dR4dY3
%     dFint4dY3=zeros(neq,neq);
%     dFint4dY3(:,1:6)=zeros(neq,6);
%     DFDUcoup43=dfdu(delta_t,struc_param,x,penal,4,3,T,T_hist,iter);
%     dFint4dY3(:,7:12)=DFDUcoup43(:,7:12);
%     dFint4dY3(:,13:14)=zeros(neq,2);
%     dFint4dY3(:,15:16)=DFDUcoup43(:,15:16);
%     dFint4dY3(:,17:18)=zeros(neq,2);
%     
%     dFext4dY3=zeros(neq,neq);
%     
%     dR4dY3=dFint4dY3-dFext4dY3;
%     
%     %dR3dY3
%     dFint3dY3=zeros(neq,neq);
%     dFint3dY3(:,1:6)=zeros(neq,6);
%     dFint3dY3(:,7:12)=KTAN(:,7:12,iter-1);
%     dFint3dY3(:,13:14)=zeros(neq,2);
%     dFint3dY3(:,15:16)=KTAN(:,15:16,iter-1);
%     dFint3dY3(:,17:18)=zeros(neq,2);
%     
%     dFext3dY3=zeros(neq,neq);
%     dFext3dY3(1:6,1:6)=eye(6);
%     dFext3dY3(13,13)=1;
%     dFext3dY3(14,14)=1;
%     dFext3dY3(17,17)=1;
%     dFext3dY3(18,18)=1;
%    
%     dR3dY3= dFint3dY3-dFext3dY3;
%     lambda(:,iter-1)=-[lambda(:,iter)'*dR4dY3]/[dR3dY3];
% %-------------------------------------------------------------------------
%    %dR4dY2
%     dFint4dY2=zeros(neq,neq);
%     dFint4dY2(:,1:6)=zeros(neq,6);
%     DFDUcoup42=dfdu(delta_t,struc_param,x,penal,4,2,T,T_hist,iter);
%     dFint4dY2(:,7:12)=DFDUcoup42(:,7:12);
%     dFint4dY2(:,13:14)=zeros(neq,2);
%     dFint4dY2(:,15:16)=DFDUcoup42(:,15:16);
%     dFint4dY2(:,17:18)=zeros(neq,2);
%     
%     dFext4dY2=zeros(neq,neq);
%     
%     dR4dY2=dFint4dY2-dFext4dY2;
%   
%    %dR3dY2
%     dFint3dY2=zeros(neq,neq);
%     dFint3dY2(:,1:6)=zeros(neq,6);
%     DFDUcoup32=dfdu(delta_t,struc_param,x,penal,3,2,T,T_hist,iter);
%     dFint3dY2(:,7:12)=DFDUcoup32(:,7:12);
%     dFint3dY2(:,13:14)=zeros(neq,2);
%     dFint3dY2(:,15:16)=DFDUcoup32(:,15:16);
%     dFint3dY2(:,17:18)=zeros(neq,2);
%     
%     dFext3dY2=zeros(neq,neq);
%     
%     dR3dY2=dFint3dY2-dFext3dY2;
%     
%      
%   %dR2Y2  
%     dFint2dY2=zeros(neq,neq);
%     dFint2dY2(:,1:6)=zeros(neq,6);
%     dFint2dY2(:,7:12)=KTAN(:,7:12,iter-2);
%     dFint2dY2(:,13:14)=zeros(neq,2);
%     dFint2dY2(:,15:16)=KTAN(:,15:16,iter-2);
%     dFint2dY2(:,17:18)=zeros(neq,2);
%     
%     dFext2dY2=zeros(neq,neq);
%     dFext2dY2(1:6,1:6)=eye(6);
%     dFext2dY2(13,13)=1;
%     dFext2dY2(14,14)=1;
%     dFext2dY2(17,17)=1;
%     dFext2dY2(18,18)=1;
%    
%     dR2dY2= dFint2dY2-dFext2dY2;
%     lambda(:,iter-2)=-[lambda(:,iter-1)'*dR3dY2+lambda(:,iter)'*dR4dY2]/[dR2dY2];
% 
% 
% %------------------------------------------------------------------------
% %dR4dY1
% dFint4dY1=zeros(neq,neq);
% dFint4dY1(:,1:6)=zeros(neq,6);
% DFDUcoup41=dfdu(delta_t,struc_param,x,penal,4,1,T,T_hist,iter);
% dFint4dY1(:,7:18)=DFDUcoup41(:,7:18);
% 
% dFext4dY1=zeros(neq,neq);
% dR4dY1=dFint4dY1-dFext4dY1;
% %dR3dY1
% dFint3dY1=zeros(neq,neq);
% dFint3dY1(:,1:6)=zeros(neq,6);
% DFDUcoup31=dfdu(delta_t,struc_param,x,penal,3,1,T,T_hist,iter);
% dFint3dY1(:,7:18)=DFDUcoup31(:,7:18);
% 
% dFext3dY1=zeros(neq,neq);
% dR3dY1=dFint3dY1-dFext3dY1;
% 
% %dR2dY1
% dFint2dY1=zeros(neq,neq);
% dFint2dY1(:,1:6)=zeros(neq,6);
% DFDUcoup21=dfdu(delta_t,struc_param,x,penal,2,1,T,T_hist,iter);
% dFint2dY1(:,7:18)=DFDUcoup21(:,7:18);
% 
% dFext2dY1=zeros(neq,neq);
% dR2dY1=dFint2dY1-dFext2dY1;
% 
% %dR1dY1
% dFint1dY1=zeros(neq,neq);
% dFint1dY1(:,7:18)=KTAN(:,7:18,iter-3);
% dFext1dY1=zeros(neq,neq);
% dFext1dY1(1:6,1:6)=eye(6);
% dR1dY1=dFint1dY1-dFext1dY1;
% 
% lambda(:,iter-3)=-[lambda(:,iter-1)'*dR3dY1+lambda(:,iter)'*dR4dY1+lambda(:,iter-2)'*dR2dY1]/[dR1dY1];

%tic()
    %DFDUcoup32=dfdu(delta_t,struc_param,x,penal,3,2,T,T_hist,iter,iter_heat);
 % DFDUcoup32_f= DFDUcoup32(freedofs,freedofs);
  % nnz(DFDUcoup32_f)
    %DFDUcoup(freedofs,freedofs)
 % toc()  
       for i=iter-1:-1:1
        RHS=zeros(1,neq);
        for k=i+1:iter
            
            DFDUcoup=dfdu(delta_t,struc_param,x,penal,k,i,T,T_hist,iter,iter_heat);

            RHS(1,freedofs)=RHS(1,freedofs)+[lambda(freedofs,k)'*DFDUcoup(freedofs,freedofs)];
        end
        lambda(freedofs,i)=-RHS(1,freedofs)/KTAN(freedofs,freedofs,i);
    end
    
    %%%%%SENSITIVITY CALCULATION THROUGH ADJOINT FORMULATION
    for i=1:iter
        %         dF_dx=dFintdx(:,:,i);
        %         lambda_i=lambda(:,i);
        dfdrho=dfdrho+(lambda(:,i)'*dFintdx(:,:,i))';
        %         dfdrho=dfdrho+(lambda_i'*dF_dx)';
        
    end
     dtestdrho=zeros(nel,1);
     NDOF=2;
     for i=1:iter
        for e=1:nel
             for I=1:4
            II=(I-1)*NDOF+1;
            IDOF(II:II+1)=(LE(e,I)-1)*NDOF+1:(LE(e,I)-1)*NDOF+2;
             end
             if e==8
                 if i==8
             lambda(IDOF,i);
             dFintdx(IDOF,e,i);
                 end
             end
             
        dtestdrho(e,1)=dtestdrho(e,1)+(lambda(IDOF,i)'*dFintdx(IDOF,e,i))';
        end     
        
    end
    
    
end
end