function dofs =  ApplyBC(nodes, condition)
dofs_x = nodes.*2-1;
dofs_y = nodes.*2;

    if condition =="roller x"
       dofs = dofs_y;
    end

    if condition =="roller y"
       dofs = dofs_x;
    end

    if condition =="fixed x"
       dofs = dofs_x;
    end

    if condition =="fixed y"
       dofs = dofs_y;
    end
end