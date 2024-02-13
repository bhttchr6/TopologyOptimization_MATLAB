function PlotBVP(XYZ,LE, BC_dofs, Fnodes, Fdir, U, flag)



nodes = size(XYZ,1);

%% extract coordinates which are fixed
BC_Nodes = unique(ceil(BC_dofs./2));
XYZ_BC = XYZ(BC_Nodes,:);

%% extract force nodes coordinates
XYZ_F = XYZ(Fnodes,:);
%if abs(Fdir) == 1
%   arrow_start = XYZ
%end

DELTA=U;
counter1=0;
for i=1:size(XYZ,1)
    for j = 1:size(XYZ,2)
        counter1=counter1+1;
        CDISP(counter1,:)=XYZ(i,j)+DELTA(counter1,1);
        
    end
end


nel = size(LE,1);
for  el=1:nel
val=(LE(el,:));
n1=val(:,1);
        n2=val(:,2);
        n3=val(:,3);
        n4=val(:,4);
        %         NEW =CDISP([2*n1+1;2*n1+2;2*n2+1;2*n2+2; 2*n3+1;2*n3+2;2*n4+1;2*n4+2],1);
        NEW =CDISP([2*n1-1;2*n1;2*n2-1;2*n2; 2*n3-1;2*n3;2*n4-1;2*n4],1);
        %         n1=val(:,4);
        %         n2=val(:,3);
        %         NEW =CDISP([2*n1+1;2*n1+2;2*n2+1;2*n2+2; 2*n2-1;2*n2;2*n1-1;2*n1],1);
        xx(:,el) = [NEW(1,1) NEW(3,1) NEW(5,1) NEW(7,1)]';
        yy(:,el) = [NEW(2,1) NEW(4,1) NEW(6,1) NEW(8,1)]';
        xx_mesh(:,el) = [NEW(1,1) NEW(3,1) NEW(5,1) NEW(7,1) NEW(1,1)]';
        yy_mesh(:,el) = [NEW(2,1) NEW(4,1) NEW(6,1) NEW(8,1) NEW(2,1)]';   
end

plot(xx_mesh,yy_mesh,'r');
hold on
if flag ==2
    plot(xx_mesh,yy_mesh,'b');
end


plot(XYZ_BC(:,1), XYZ_BC(:,2), 'o', 'MarkerSize', 10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','r');
hold on
plot(XYZ_F(:,1), XYZ_F(:,2), '^', 'MarkerSize', 10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');

set(gcf, 'Renderer', 'painters');
set(gcf,'Position',[100,100,500,500])
xlabel('X(mm)', 'FontSize',14,'FontWeight','bold','Color','k')
ylabel('Y(mm)', 'FontSize',14,'FontWeight','bold','Color','k')
box off
axis off
axis equal
set(gcf,'color','w');




%legend('Design Domain', 'Location', 'northeast');


end