function [cauchy,cauchy_d1, D]= LinearElasticity(Ys, dYsdrho,F)

size_F = size(F,2);

%% material properties
%E = 200*10e3; % 200*10^3 N/mm^2
pr = 0.3;
mu = Ys/(2*(1 + pr));
lambda = Ys*pr/((1 + pr)*(1 - 2*pr));

%% derivative of mu and lambda
dmudrho = dYsdrho/(2*(1 + pr));
dlambdadrho = dYsdrho*pr/((1 + pr)*(1 - 2*pr));

%% Identity matrix
I = eye(size_F);

%% Green Strain
E = 0.5*(F'*F-I);
%% the 1st Piola Kirchoff Stress
%P = mu*(F + F' - 2*I) + lambda*trace(F - I)*I;
P = F*(2*mu*E + lambda*trace(E)*I);
%% cauchy stress
J = det(F);
cauchy = (1/J)*P*F';

%% derivative of cauchy
dPdrho = F*(2*dmudrho*E + dlambdadrho*trace(E)*I);
cauchy_d1 = (1/J)*dPdrho*F';

%% Stiffness matrix
D = Ys/((1+pr)*(1-2*pr))*[1-pr pr   0;
                          pr  1-pr 0;
                          0   0    0.5-pr];

end