function [XYZ, LE, DOMPARAM] = Domain(nelx, nely, L, H, ndof)
n_nodes_x = nelx+1;
n_nodes_y = nely+1;

%% aspect ratio
a = L/nelx;
b = H/nely;
%nodal coordinates
N=0;
for j=1:n_nodes_x
    for i=1:n_nodes_y
        N=N+1;
        XYZ(N,:)=[(j-1)*a (nely-(i-1))*b];

    end
end

%element connectivity
counter_le=0;
for elx = 1:nelx
    for ely = 1:nely
        counter_le=counter_le+1;
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx +ely;
        n3=(nely+1)*(elx-1)+ely+1;
        n4=(nely+1)* elx +ely+1;
        LE(counter_le,:)=[n3 n4 n2 n1];
    end
end

%%% number of equations
nnodes = (nelx+1)*(nely+1);
neq = nnodes*ndof;

DOMPARAM.nnodes = nnodes;
DOMPARAM.neq = neq;
DOMPARAM.nel = size(LE,1);
DOMPARAM.eps = a;
DOMPARAM.nelx = nelx;
DOMPARAM.nely = nely;
DOMPARAM.L = L;
DOMPARAM.H = H;
end
