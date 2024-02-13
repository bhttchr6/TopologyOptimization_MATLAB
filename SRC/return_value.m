% FUNCTION TO EXTRACT VALUE AT A PARTICULAR TIME

function [m_U,m_theta,dL_dx]=return_value(nelx,nely,time_spec,delta_t,x,penal,sens)

[disp,force_val,temp,time,dL_dx]=cyclic_loading(nelx,nely,time_spec,delta_t,x,penal,sens);
data_val=[disp;force_val;temp;time];
for i=1:size(data_val,2)
    if data_val(4,i)<=time_spec
        counter=i;
    end
end
m_U=data_val(1,counter);
m_theta=data_val(2,counter);
end