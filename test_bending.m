%discrete_shell_func = [];
%[V,F] = subdivided_sphere(3);
%[V,F] = create_irregular_grid(50,50);
%[V,F] = create_regular_grid(10,10);
[V,F] = load_mesh('~/Dropbox/models/cartoon-elephant.ply');

%F = [1 2 3;2 1 4];
%V = [0 -1;0 1;-1 0;1 0];

V(:,end+1:3) = 0;
[EF,EI,E,~,EJ] = edge_flaps(F);
keep = all(EF>0,2);
EF = EF(keep,:);
EI = EI(keep,:);
E = E(keep,:);
EJ = EJ(keep,:);
EK = F(sub2ind(size(F),EF,EI));
% Scharnier
S = [E EK];

V0 = V;
[th0,w0] = dihedral_angle(V0,S);
if isempty(discrete_shell_func)
  tic;
  discrete_shell_func = per_element_energy(@hinge,V,S,'Name','discrete_shell','Constants',{th0,w0});%,'UseFile',false);
  toc
end

tic;
[f,G,H] = discrete_shell_func(V,S,th0,w0);
toc
[EV,ED] = eigs(H+H',repdiag(massmatrix(V,F),size(V,2)),15,'sm');
% huh... nullspace of bending on a sphere appears to be 7-dimensional

for e = 8:size(EV,2)
for t = interp1([0 1 2],[0 1 0],linspace(0,2));
  tsurf(F,V+t*reshape(EV(:,e),size(V)),'CData',t*normrow(reshape(EV(:,e),size(V))),fphong);
  caxis([0 max(normrow(reshape(EV(:,e),size(V))))]);
  axis equal;
  drawnow;
end
end


function [th,w] = dihedral_angle(V,S)
  % V  4 by 3 list of vertex positions, triangle then flaps
  N1 = cross(V(S(:,1),:)-V(S(:,3),:),V(S(:,2),:)-V(S(:,3),:),2);
  N2 = cross(V(S(:,2),:)-V(S(:,4),:),V(S(:,1),:)-V(S(:,4),:),2);
  % Twice the area
  A1 = sqrt(sum(N1.^2,2));
  A2 = sqrt(sum(N2.^2,2));
  %N1 = N1./A1;
  %N2 = N2./A2;
  % Unit edge vector
  EV = V(S(:,2),:)-V(S(:,1),:);
  l2 = sum(EV.^2,2);
  EV = EV./sqrt(l2);
  th = pi-atan2(dot(cross(N1,N2,2),EV,2),dot(N1,N2,2));
  w = 3*l2./(A1+A2);
end

function f = hinge_nodal(V,V0)
  % V  4 by 3 list of vertex positions, triangle then flaps
  % V0  4 by 3 same vertices at rest
  th = dihedral_angle(V,[1 2 3 4]);
  [th0,w0] = dihedral_angle(V0,[1 2 3 4]);
  f = w0*(th-th0).^2;
end
function f = hinge(V,th0,w0)
  % V  4 by 3 list of vertex positions, triangle then flaps
  % V0  4 by 3 same vertices at rest
  th = dihedral_angle(V,[1 2 3 4]);
  f = w0*(th-th0).^2;
end
