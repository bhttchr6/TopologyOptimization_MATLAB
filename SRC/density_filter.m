function [ W] = density_filter( XYZ,LE,rmin,nel)

for e=1:nel
    elXY=XYZ(LE(e,:),:); % elemental nodal coordinates
    X=elXY(:,1); % elemnt X coordinate
    Y=elXY(:,2); % element Y coordinate
    x_c(e)=mean(X);
    y_c(e)=mean(Y);
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



        %           rmin=35*0.0005*x_c1;
        w(element_1,element_2)=max(0,rmin-R);


    end

    wnorm=sum(w(element_1,:));
    W(element_1,:) = w(element_1,:)/wnorm;


end


end
