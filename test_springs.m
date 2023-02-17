
n = 20;
X = interp1([0;1],[0 0 1;1 0 1],linspace(0,1,n));
E = [1:size(X,1)-1;2:size(X,1)]';

k = repmat(1e8,size(E,1),1);
R = sqrt(sum((X(E(:,2),:)-X(E(:,1),:)).^2,2));
%
%
x = X;
spring_func = per_element_energy(@spring_3d,x,E,'Name','spring_3d','Constants',{k,R},'Nodal',[false false]);
%[f,G,H] = spring_func(X*2,E,k,R);
%func_nodal = per_element_energy(@spring_3d_nodal,x,E,'Name','spring_3d_nodal','Constants',{k,X},'Nodal',[false true]);
%[f,G,H] = func_nodal(X*2,E,k,X)

x = X;
x0 = x;
x1 = x;
b = 1;
bc = X(b,:);

M = repdiag(sparse(E,E,[R R]/2,size(X,1),size(X,1)),size(X,2));
h = 1/300;

g = zeros(size(X));
g(:,3) = -9.8;

for outer = 1:1000
  x1 = x0;
  x0 = x;

% only depends on x,x0,x1
f = @(x) spring_func(x,E,k,R) +  ...
  0.5 * 1/8/h^2*(x(:)- 2*x0(:) + x1(:))'*M*(x(:) - 2*x0(:) + x1(:)) + ...
  -x(:)'*M*g(:);
for newton_iter = 1:100
  [f_sp,G_sp,H_sp] = spring_func(x,E,k,R);
  G_mo =  1/8*M*(x(:)-2*x0(:)+x1(:));
  H_mo =  1/8*M;
  H_ext =  h^2*sparse(0);
  G_ext = -h^2*M*g(:);
  G =  h^2*G_sp(:) + G_mo + G_ext;
  H =  h^2*H_sp + H_mo + H_ext;

  %dx = reshape(-G,size(x));
  dx = reshape(min_quad_with_fixed( ...
    0.5*H, ...
    G, ...
    b+[0:size(x,2)-1]*size(x,1),bc(:)*0),size(x));
  if norm(dx)<1e-10
    break;
  end
  [t,x,fx] = backtracking_line_search(f,x,reshape(G,size(x)),dx,0.3,0.5,30);
    %plot_edges(x,E,'-ok','LineWidth',0.5)
    %axis equal;
    %axis([-1 1 -1 1 -1 1]);
    %drawnow;
  if t == 0
    break;
  end
  if newton_iter == 100
    warning('newton did not converge');
  end
end
plot_edges(x,E,'-ok','LineWidth',1)
axis equal;
axis([-1 1 -1 1 -1 1]);
drawnow;

end



function f = spring_3d(x,k,R)
  % Inputs:
  %   x  2 by 3 list of vertex positions
  %   k  scalar spring coefficient
  %   R  scalar spring rest length
  %
  r = sqrt(sum((x(2,:)-x(1,:)).^2));
  f = k*(r-R)^2;
end

function f = spring_3d_nodal(x,k,X)
  % Inputs:
  %   x  2 by 3 list of vertex positions
  %   k  scalar spring coefficient
  %   X  2 by 3 list of vertex positions
  %
  r = sqrt(sum((x(2,:)-x(1,:)).^2));
  R = sqrt(sum((X(2,:)-X(1,:)).^2));
  f = k*(r-R)^2;
end
