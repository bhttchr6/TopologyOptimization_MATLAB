function [dNdx,B, det] = SHAPEL2DTRI(C)

% C is the current coordinates [x1 y1;x2 y2;x3 y3];

x2=C(2,1);
x1=C(1,1);
x3=C(3,1);
y1=C(1,2);
y2=C(2,2);
y3=C(3,2);

det=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
B=(1/det)*[y2-y3  0    y3-y1   0    y1-y2     0;
    0   x3-x2   0   x1-x3    0    x2-x1;
    x3-x2 y2-y3 x1-x3 y3-y1  x2-x1   y1-y2];


j_inv=(1/det)*[  y3-y1    -(y2-y1);
    -(x3-x1)      x2-x1];

dNdx=j_inv*[-1 1 0;-1 0 1];

end