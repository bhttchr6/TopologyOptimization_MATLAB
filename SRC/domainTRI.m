function domainTRI(XYZ,LE,rhoVec,U,mes,den,outeriter,kktnorm,vol,U_end,penal)
global num BC
tri=1;
% number of elements
nel=size(LE,1);

DELTA=U;
counter1=0;
for i=1:size(XYZ,1)
    for j = 1:size(XYZ,2)
        counter1=counter1+1;
        CDISP(counter1,:)=XYZ(i,j)+DELTA(counter1,1);

    end
end

if tri==1;


    parfor el = 1:nel
        el;
        val=(LE(el,:));
        n1=val(:,1);
        n2=val(:,2);
        n3=val(:,3);
        NEW =CDISP([2*n1-1;2*n1;2*n2-1;2*n2;2*n3-1;2*n3],1);

        xx(:,el) = [NEW(1,1) NEW(3,1) NEW(5,1)]';
        xx_mesh(:,el)=[NEW(1,1) NEW(3,1) NEW(5,1) NEW(1,1)]'
        yy(:,el)= [NEW(2,1) NEW(4,1) NEW(6,1)]';
        yy_mesh(:,el)= [NEW(2,1) NEW(4,1) NEW(6,1) NEW(2,1)]';
    end
end





if tri==0

    parfor el=1:nel
        val=(LE(el,:));
        n1=val(:,4);
        n2=val(:,3);
        NEW =CDISP([2*n1+1;2*n1+2;2*n2+1;2*n2+2; 2*n2-1;2*n2;2*n1-1;2*n1],1);

        xx(:,el) = [NEW(1,1) NEW(3,1) NEW(5,1) NEW(7,1)]';
        yy(:,el) = [NEW(2,1) NEW(4,1) NEW(6,1) NEW(8,1)]';
        xx_mesh(:,el) = [NEW(1,1) NEW(3,1) NEW(5,1) NEW(7,1) NEW(1,1)]';
        yy_mesh(:,el) = [NEW(2,1) NEW(4,1) NEW(6,1) NEW(8,1) NEW(2,1)]';
    end
end
fixed_nodes=BC';
% fixed_nodes=[24:68 331:372];
% fixed_nodes=[9:32 53 82:93  518 519 554:575 462:496];
% fixed_nodes=[9 82 83 84 85 86 87 88 89 90 91 92 93 554 555 556 557 558 559 560 ...
%       561 562 563 564 565 566 567 568 569 570 571 ...
%       572 573 574 575 10 11 12 13 14 15 16 17 18 19 20 ...
%       21 22 23 24 25 26 27 28 29 30 31 32 ...
%       462 463 464 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 480 481 ...
%       482 483 484 485 486 487 488 489 490 ...
%       491 492 493 494 495 496 ];
% fixed_nodes=[ 53 82 89 90 91 92 93  518 519 566 567 568 569 570 571 ...
%       572 573 574 575  17 18 19 20 ...
%       21 22 23 24 25 26 27 28 29 30 31 32 ...
%       471 472 473 474 475 476 477 478 479 480 481 ...
%       482 483 484 485 486 487 488 489 490 ...
%       491 492 493 494 495 496];
fixed_nodes=sort(fixed_nodes);
for i=1:size(fixed_nodes,2)
    XYZ_P(i,:)=(XYZ(fixed_nodes(i),:));
end
% Mesh Plotting
if mes==1
    figure(2)
    plot(xx_mesh,yy_mesh,'r')

    hold on
    plot (XYZ_P(:,1),XYZ_P(:,2),'*')
    % text(XYZ(462,1)+550,XYZ(462,2)+550,'\bf Fixed Boundary ')
    % text(XYZ(474,1)-70,XYZ(474,2)-70,'\bf Fixed Boundary ')
    % hold on
    % % plot(XYZ(1833,1),XYZ(1833,2),'*')
    %
    % text(XYZ(1833,1)-5,XYZ(1833,2),' \bf Force\bullet','HorizontalAlignment','right')
    % hold on
    % % plot(XYZ(52,1),XYZ(52,2),'*')
    %
    % text(XYZ(52,1)-25,XYZ(52,2)+10,' \bullet \bf Deflection')
    colormap("jet");
    set(gcf,'color','w');
    grid off
    axis off
    axis equal
end

% Density distribution plotting
if den==1
    h1=figure(2)
    patch(xx,yy,-rhoVec)
    set(gcf,'color','w');
    shading flat
    colormap("jet");
    grid off
    axis off
    %title([num2str(outeriter) ' ' 'kktnorm:' num2str(kktnorm) ' ''volfrac:' num2str(vol) ' ' num2str(U_end) ' ' num2str(penal) ' ' num2str(num)]);
    axis equal
    drawnow;
    saveas(h1,sprintf('FIGURE%d.png',outeriter))
end
end