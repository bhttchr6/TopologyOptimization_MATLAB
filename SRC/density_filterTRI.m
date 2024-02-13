function [ W] = density_filterTRI( XYZ,LE,rmin,nel)


% nel=nelx*nely;
% nnodes=(nelx+1)*(nely+1);
% neq=nnodes*2;
% n_nodes_x=nelx+1;
% n_nodes_y=nely+1;


% coordinates of each element
%  N=0;
% for j=1:n_nodes_x
%     for i=1:n_nodes_y
%         N=N+1;
%         XYZ(N,:)=[(j-1)*a (nely-(i-1))*b];
%
%     end
% end


% centroid of each element
% for elx=1:nelx
%     for ely=1:nely
%
%
%  n1 = (nely+1)*(elx-1)+ely;
%  n2 = (nely+1)* elx +ely;
%  n3=(nely+1)*(elx-1)+ely+1;
%  n4=(nely+1)* elx +ely+1;
%       x1=XYZ(n1,:);
%       x2=XYZ(n2,:);
%       x3=XYZ(n3,:);
%       x4=XYZ(n4,:);
%       x_c(ely,elx)=((x2(1)-x1(1))/2)+(elx-1)*a;
%       y_c(ely,elx)=((x2(2)-x4(2))/2)+(ely-1)*b;
%
%
%
%     end
% end
for e=1:nel
    elXY=XYZ(LE(e,:),:); % elemental nodal coordinates
    X=elXY(:,1); % elemnt X coordinate
    Y=elXY(:,2); % element Y coordinate
    x_c(e)=mean(X);
    y_c(e)=mean(Y);
    %     rmin(e)=sqrt((X(1)-X(2))^2+((Y(1)-Y(2))^2));
end

element_1=0;

for el_1=1:nel

    element_1=element_1+1;

    x_c1=x_c(el_1);
    y_c1=y_c(el_1);
    element_2=0;
    for el_2=1:nel

        element_2=element_2+1;
        x_c2=x_c(el_2);
        y_c2=y_c(el_2);
        R=sqrt((y_c2-y_c1)^2+(x_c2-x_c1)^2);
        w(element_1,element_2)=max(0,rmin-R);

    end

    wnorm=sum(w(element_1,:));
    W(element_1,:) = w(element_1,:)/wnorm;


end




end
