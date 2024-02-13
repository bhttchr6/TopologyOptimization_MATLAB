function force = ForceEval(FCond, GDof)

%% nodes at which forces are applied and the corresponding angles


force = zeros(GDof,1);

Force_mat = diag(FCond(:,3))*[cosd(FCond(:,4)) sind(FCond(:,4))];
ForceMagDir = reshape(Force_mat.',1,[])';
Fdofs = [FCond(:,1).*2-1 FCond(:,1).*2];
force(Fdofs,1) = ForceMagDir;

end