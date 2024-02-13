function Force =  ApplyForce(FCond, DOMPARAM,  condition)
neq = DOMPARAM.neq;
nodes = FCond.ForceNodes;
Magnitude = FCond.ForceMagnitude;
Force = zeros(neq,1);
dofs_x = nodes.*2-1;
dofs_y = nodes.*2;

   
    if condition =="along x"
       dofs = dofs_x;
    end

    if condition =="along y"
       dofs = dofs_y;
    end

    Force(dofs,1) = Magnitude;
end